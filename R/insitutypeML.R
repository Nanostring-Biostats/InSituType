#' Classify cells based on reference profiles
#' 
#' Supervised classification of cells. Each cell is assigned to the cell type 
#'  under which its observed expression profile is most likely. 
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell. Can be provided 
#' @param bg Expected background
#' @param cohort Vector of cells' cohort memberships
#' @param reference_profiles Matrix of expression profiles of pre-defined clusters,
#'  e.g. from previous scRNA-seq. These profiles will not be updated by the EM algorithm.
#'  Colnames must all be included in the init_clust variable.
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param align_genes Logical, for whether to align the counts matrix and the reference_profiles by gene ID.
#' @return A list, with the following elements:
#' \enumerate{
#' \item clust: a vector given cells' cluster assignments
#' \item prob: a vector giving the confidence in each cell's cluster
#' \item profiles: Matrix of clusters' mean background-subtracted profiles
#' \item logliks: Matrix of cells' log-likelihoods under each cluster. Cells in rows, clusters in columns.
#' }
#' @export
insitutypeML <- function(counts, neg = NULL, bg = NULL, cohort = NULL, reference_profiles, nb_size = 10, align_genes = TRUE) {
  
  # get vector of expected background:
  if (is.null(bg) & is.null(neg)) {
    stop("Must provide either bg or neg")
  }
  # infer bg from neg if needed
  if (is.null(bg) & !is.null(neg)) {
      s <- Matrix::rowMeans(counts)
      bgmod <- stats::lm(neg ~ s - 1)
      bg <- bgmod$fitted
      names(bg) = rownames(counts)
  }
  # accept a single value of bg if input by user:
  if (length(bg) == 1) {
    bg <- rep(bg, nrow(counts))
    names(bg) <- rownames(counts)
  }
  
  # align genes:
  if (align_genes) {
    sharedgenes <- intersect(rownames(reference_profiles), colnames(counts))
    lostgenes <- setdiff(colnames(counts), rownames(reference_profiles))
    
    # subset:
    counts <- counts[, sharedgenes]
    reference_profiles <- reference_profiles[sharedgenes, ]
    
    # warn about genes being lost:
    if ((length(lostgenes) > 0) & length(lostgenes) < 50) {
      message(paste0("The following genes in the count data are missing from reference_profiles and will be omitted from cell typing: ",
                     paste0(lostgenes, collapse = ",")))
    }
    if (length(lostgenes) > 50) {
      message(paste0(length(lostgenes), " genes in the count data are missing from reference_profiles and will be omitted from cell typing"))
    }
  }
  
  # prep cohort vector:
  if (is.null(cohort)) {
    cohort <- rep("all", length(bg))
  }
  
  # get logliks
  logliks <- apply(reference_profiles, 2, function(ref) {
    lldist(x = ref,
           mat = counts,
           bg = bg,
           size = nb_size)
  })
  
  # update logliks based on frequencies within cohorts:
  logliks <- update_logliks_with_cohort_freqs(logliks = logliks, 
                                              cohort = cohort, 
                                              minfreq = 1e-4, 
                                              nbaselinecells = 100) 
  
  # get remaining outputs
  clust <- colnames(logliks)[apply(logliks, 1, which.max)]
  names(clust) <- rownames(logliks)
  
  probs <- logliks2probs(logliks)
  prob <- apply(probs, 1, max)
  names(prob) <- names(clust)
  profiles <- Estep(counts, 
                    clust = clust,
                    neg = neg)
  # aligns profiles and logliks, removing lost clusters:
  logliks_from_lost_celltypes <- logliks[, !is.element(colnames(logliks), unique(clust)), drop = FALSE]
  logliks <- logliks[, is.element(colnames(logliks), clust), drop = FALSE]
  profiles <- profiles[, colnames(logliks), drop = FALSE]
  
  out = list(clust = clust,
             prob = prob,
             profiles = profiles,
             logliks = round(logliks, 4),
             logliks_from_lost_celltypes = round(logliks_from_lost_celltypes, 4))
  return(out)    
  break()
}