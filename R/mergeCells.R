#' Merge cell types in a clustering result
#'
#' Take a user-defined list of cells types to rename/combine, then re-compute
#'  cluster assignments and probabilities under the merged cell types.
#' @param merges A named vector in which the elements give new cluster names and
#'  the names give old cluster names. OK to omit cell types that aren't being merged.
#' @param to_delete Cluster names to delete. All cells assigned to these clusters 
#'  will be reassigned to the next best cluster. 
#' @param probs Matrix of probabilities output by insitutype
#' @param profiles Profiles matrix output by insitutype
#' @return A list with two elements:
#' \enumerate{
#' \item clust: a vector of cluster assignments
#' \item probs: a matrix of probabilities of all cells (rows) belonging to all clusters (columns)
#' }
#' @export
#' @examples
#' # define a "merges" input:
#' merges = c("macrophages" = "myeloid", "monocytes" = "myeloid", "mDC" = "myeloid",
#'              "B-cells" = "lymphoid")
#' # define clusters to delete:
#' to_delete =  c("a", "f") 
mergeCells <- function(merges = NULL, to_delete = NULL, probs, profiles) {
  
  # input checks:
  #if (is.null(merges) & is.null(to_delete)) {
  #  stop("Must provide a value foe either merges or to_delete")
  #}
  
  # check that merges and to_delete names are all in logliks:
  if (any(!is.element(names(merges), colnames(probs)))) {
    mismatch <- setdiff(names(merges), colnames(probs))
    stop(paste0("The following user-provided cluster name(s) in the merges argument are missing from colnames(probs): ",
                paste0(mismatch, collapse = ", ")))
  }
  if (any(!is.element(to_delete, colnames(probs)))) {
    mismatch <- setdiff(to_delete, colnames(probs))
    stop(paste0("The following user-provided cluster name(s) in the to_delete argument are missing from colnames(probs): ",
                paste0(mismatch, collapse = ", ")))
  }
  if (length(setdiff(colnames(probs), to_delete)) == 0) {
    stop("The to_delete argument is asking for all clusters to be deleted.")
  }
  
  # convert probs to logliks:
  logliks <- probs2logliks(probs) 
    
  # delete those called for:
  logliks <- logliks[, !is.element(colnames(logliks), to_delete)]
  
  # get logliks under merged categories: each cell's "new" loglik in a merged cell type is
  #  its best loglik under the "old" celltype.
  newlogliks <- matrix(NA, nrow(logliks), length(unique(merges)),
                       dimnames = list(rownames(logliks), unique(merges)))
  newlogliks <- sapply(unique(merges), function(newname) {
    oldnames <- names(merges)[merges == newname]
    newlogliks[, newname] = apply(logliks[, oldnames, drop = FALSE], 1, max)
  })
  if (length(newlogliks) > 0) {
    newlogliks <- cbind(newlogliks, logliks[, setdiff(colnames(logliks), names(merges)), drop = FALSE])
  } else {
    newlogliks <- logliks
  }

  ## convert to probs:
  probs <- logliks2probs(newlogliks)
  # flag cells with 0 probability in any of the remaining cell types:
  lostcells <- is.na(probs[, 1])
  
  # get new cluster assignments:
  clust <- rep(NA, nrow(probs))
  names(clust) <- rownames(probs)
  clust[!lostcells] <- colnames(probs)[apply(probs[!lostcells, , drop = FALSE], 1, which.max)]
  out <- list(clust = clust, probs = round(probs, 3))  # (rounding probs to save memory)
  return(out)
}


# get a logliks matrix from a probabilities matrix
probs2logliks <- function(probs) {
  return(log(probs))
}


#' convert logliks to probabilities
#' 
#' From cell x cluster log-likelihoods, calculate cell x cluster probabilities
#' @param logliks Matrix of loglikelihoods, as output by insitytupe. Cells in rows, clusters in columns.
#' @return A matrix of probabilities, in the same dimensions as logliks. 
#' @export 
logliks2probs <- function(logliks) {
  templogliks <- sweep(logliks, 1, apply(logliks, 1, max ), "-" )
  # get on likelihood scale:
  liks <- exp(templogliks)
  # convert to probs
  probs <- sweep(liks, 1, rowSums(liks), "/")
  return(probs)
}


