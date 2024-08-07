#' Merge cell types in a clustering result
#'
#' Take a user-defined list of cells types to rename/combine, then re-compute
#' cluster assignments and probabilities under the merged cell types.
#' @param assay_type Assay type of RNA, protein (default = "rna")
#' @param merges A named vector in which the elements give new cluster names and
#'   the names give old cluster names. OK to omit cell types that aren't being
#'   merged.
#' @param to_delete A vector of cluster names to delete. All cells assigned to
#'   these clusters will be reassigned to the next best cluster.
#' @param subcluster A list, where each element's name is a cell type to
#'   subcluster, and the element itself is the cluster number(s) to use. E.g.
#'   list("macrophages" = 2, "cancer" = 2:3)
#' @param logliks Matrix of log-likelihoods output by insitutype, cells in rows,
#'   clusters in columns
#' @param counts Counts matrix, cells * genes. Only needed if subclustering is
#'   run.
#' @param neg Vector of mean negprobe counts per cell. Only needed if
#'   subclustering is run.
#' @param bg Expected background. Optional, and only used if subclustering is
#'   run.
#' @param cohort Vector of cells' cohort memberships. Optional, and only needed
#'   if subclustering is run.
#' @return A list with two elements: \enumerate{ \item clust: a vector of
#'   cluster assignments \item prob: Vector of posterior probabilities for each
#'   cell type \item logliks: a matrix of probabilities of all cells (rows)
#'   belonging to all clusters (columns) \item profiles: a matrix of the average
#'   background-subracted profile of each cell type after
#'   merging/deleting/subclustering }
#' @export
#' @examples
#' #example merges argument:
#' merges = c("macrophages" = "myeloid",  # merge 3 clusters
#'            "monocytes" = "myeloid",
#'            "mDC" = "myeloid",
#'            "B-cells" = "lymphoid")    # just rename 1 cluster
#' # example to_delete argument:
#' to_delete = c("neutrophils")
#' # example subcluster argument:
#' subcluster = list("Myofibroblast" = 2:3)
refineClusters <- function(assay_type = c("rna", "protein"), 
                           merges = NULL, to_delete = NULL, subcluster = NULL, 
                           logliks,
                           counts = NULL, 
                           neg = NULL, bg = NULL, 
                           cohort = NULL) {
  assay_type <- match.arg(tolower(assay_type), c("rna", "protein"))
  
  # check that provided cell names are all in logliks:
  if (any(!is.element(names(merges), colnames(logliks)))) {
    mismatch <- setdiff(names(merges), colnames(logliks))
    stop(paste0("The following user-provided cluster name(s) in the merges argument are missing from colnames(logliks): ",
                paste0(mismatch, collapse = ", ")))
  }
  if (any(!is.element(to_delete, colnames(logliks)))) {
    mismatch <- setdiff(to_delete, colnames(logliks))
    stop(paste0("The following user-provided cluster name(s) in the to_delete argument are missing from colnames(logliks): ",
                paste0(mismatch, collapse = ", ")))
  }
  if (any(!is.element(names(subcluster), colnames(logliks)))) {
    mismatch <- setdiff(names(subcluster), colnames(logliks))
    stop(paste0("The following user-provided cluster name(s) in the merges argument are missing from colnames(logliks): ",
                paste0(mismatch, collapse = ", ")))
  }
  if (length(setdiff(colnames(logliks), to_delete)) == 0) {
    stop("The to_delete argument is asking for all clusters to be deleted.")
  }
  # check that subcluster data is available:
  if (!is.null(subcluster)) {
    if (is.null(counts)) {
      stop("Must provide counts data to subcluster")
    }
    if (is.null(neg)) {
      stop("Must provide neg vector to subcluster")
    }
  }

  # delete those called for:
  logliks <- logliks[, !is.element(colnames(logliks), to_delete)]
  
  # prevent merging into existing names:
  duplicatednames <- intersect(merges, colnames(logliks))
  merges[is.element(merges, duplicatednames)] <- paste0(merges[is.element(merges, duplicatednames)], ".new")
  
  # get logliks under merged categories: each cell's "new" loglik in a merged cell type is
  #  its best loglik under the "old" celltype.
  newlogliks <- matrix(NA, nrow(logliks), length(unique(merges)),
                       dimnames = list(rownames(logliks), unique(merges)))
  newlogliks <- sapply(unique(merges), function(newname) {
    oldnames <- names(merges)[merges == newname]
    newlogliks[, newname] <- apply(logliks[, oldnames, drop = FALSE], 1, max, na.rm = TRUE)
  })
  if (length(newlogliks) > 0) {
    newlogliks <- cbind(newlogliks, logliks[, setdiff(colnames(logliks), names(merges)), drop = FALSE])
  } else {
    newlogliks <- logliks
  }
  
  # get new cluster assignments:
  clust <- colnames(newlogliks)[apply(newlogliks, 1, which.max)]
  names(clust) <- rownames(newlogliks)

  ## perform subclustering:
  # subclustering logic:
  # - run unsupervised clustering of the selected cell type 
  # - record the subcluster logliks for the selected cells
  # - for unselected cells, propagate the original supercluster loglik to the subclusters (to prevent unselected cells joining the subclusters)
  # - unselected cells keep their cell type. selected cells go to whichever subcluster gives them the greatest loglik
  for (name in names(subcluster)) {
    message(paste0("Subclustering ", name))
    use <- which(colnames(newlogliks)[apply(newlogliks, 1, which.max)] == name)
    # run insitutype on just the named cell type:
    temp <- insitutype(x = counts[use, ],
                       assay_type = assay_type,
                       neg = neg[use],
                       bg = bg[use],
                       cohort = cohort[use],
                       n_clusts = subcluster[[name]],
                       n_starts = 3, n_benchmark_cells = 5000,
                       n_phase1 = 2000, n_phase2 = 10000, n_phase3 = 20000,
                       n_chooseclusternumber = 2000)
    
    # make logliks matrix for all cells vs. the new clusters, with cells outside 
    # the selected cell type retaining their original loglik for the cluster
    subclustlogliks <- matrix(rep(newlogliks[, name], ncol(temp$logliks)), nrow(counts))
    rownames(subclustlogliks) <- rownames(counts)
    colnames(subclustlogliks) <- colnames(temp$logliks)
    # for cells with subclustering results, overwrite the old logliks:
    subclustlogliks[rownames(temp$logliks), colnames(temp$logliks)] <- temp$logliks  
    # better names:
    colnames(subclustlogliks) <- paste0(name, "_", seq_len(ncol(subclustlogliks)))
    
    # safeguard in case we've created a cell type name that already exists:
    if (any(is.element(colnames(subclustlogliks), colnames(newlogliks)))) {
      colnames(subclustlogliks) <- paste0(colnames(subclustlogliks), "subcluster")
    }
    
    # update logliks matrix:
    newlogliks <- newlogliks[, setdiff(colnames(newlogliks), name)]
    newlogliks <- cbind(newlogliks, subclustlogliks)
    
    # update clust for the subclustered cells:
    clust[use] <- colnames(subclustlogliks)[apply(subclustlogliks[use, ], 1, which.max)]
  }
  
  # get new posterior probs:
  probs <- logliks2probs(newlogliks)
  prob <- apply(probs, 1, max)
  names(prob) <- names(clust)
  
  # re-calculate profiles if available:
  profiles <- NULL
  sds <- NULL
  if (!is.null(counts) && !is.null(neg)) {
    profiles_info <- Estep(counts = counts,
                           clust = clust,
                           neg = neg,
                           assay_type=assay_type)
    profiles <- profiles_info$profiles
    sds <- profiles_info$sds
    
  }
  # aligns profiles and logliks, removing lost clusters:
  logliks_from_lost_celltypes <- newlogliks[, !is.element(colnames(newlogliks), unique(clust)), drop = FALSE]
  newlogliks <- newlogliks[, is.element(colnames(newlogliks), clust), drop = FALSE]
  profiles <- profiles[, colnames(newlogliks), drop = FALSE]
  
  if(identical(tolower(assay_type), "protein")){
    sds <- sds[, colnames(newlogliks), drop = FALSE]
  }
  
  if(identical(tolower(assay_type), "rna")){
    sds <- NULL
  }
  out <- list(clust = clust, prob = prob, logliks = round(newlogliks, 4), # (rounding logliks to save memory)
              profiles = profiles, sds=sds, logliks_from_lost_celltypes = round(logliks_from_lost_celltypes, 4))  
  return(out)
}


#' Get a logliks matrix from a probabilities matrix
#'
#' @param probs probability matrix
#'
#' @return log-transformed matrix
#' @examples 
#' a <- runif(10)
#' probs2logliks(a/sum(a))
probs2logliks <- function(probs) {
  return(log(probs))
}


#' convert logliks to probabilities
#' 
#' From cell x cluster log-likelihoods, calculate cell x cluster probabilities
#' @param logliks Matrix of loglikelihoods, as output by insitytupe. Cells in rows, clusters in columns.
#' @return A matrix of probabilities, in the same dimensions as logliks. 
#' @examples 
#' data("mini_nsclc")
#' unsup <- insitutype(
#'  x = mini_nsclc$counts,
#'  neg = Matrix::rowMeans(mini_nsclc$neg),
#'  n_clusts = 8,
#'  n_phase1 = 200,
#'  n_phase2 = 500,
#'  n_phase3 = 2000,
#'  n_starts = 1,
#'  max_iters = 5,
#'  assay_type="RNA"
#' ) # choosing inadvisably low numbers to speed the vignette; using the defaults in recommended.
#' logliks2probs(unsup$logliks)
#' 
logliks2probs <- function(logliks) {
  templogliks <- sweep(logliks, 1, apply(logliks, 1, max, na.rm = TRUE), "-")
  # get on likelihood scale:
  liks <- exp(templogliks)
  # convert to probs
  probs <- sweep(liks, 1, rowSums(liks, na.rm = TRUE), "/")
  return(probs)
}
