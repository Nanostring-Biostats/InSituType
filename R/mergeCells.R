#' Merge cell types in a clustering result
#'
#' Take a user-defined list of cells types to rename/combine, then re-compute
#'  cluster assignments and probabilities under the merged cell types.
#' @param merges A named vector in which the elements give new cluster names and
#'  the names give old cluster names. OK to omit cell types that aren't being merged.
#' @param probs Matrix of probabilities
#' @return A list with two elements:
#' \enumerate{
#' \item clust: a vector of cluster assignments
#' \item probs: a matrix of probabilities of all cells (rows) belonging to all clusters (columns)
#' }
#' @export
#' @examples
#' # define a "merges" input:
#' merges <- c("macrophages" = "myeloid", "monocytes" = "myeloid", "mDC" = "myeloid",
#'              "B-cells" = "lymphoid")
mergeCells <- function(merges, probs) {

  # convert probs to logliks:
  logliks <- probs2logliks(probs) 
    
  # check that merges names are all in logliks:
  if (any(!is.element(names(merges), colnames(logliks)))) {
    mismatch <- setdiff(names(merges), colnames(logliks))
    stop(paste0("The following user-provided cluster name(s) are missing from colnames(logliks): ",
                paste0(mismatch, collapse = ", ")))
  }

  # get logliks under merged categories: each cell's "new" loglik in a merged cell type is
  #  its best loglik under the "old" celltype.
  newlogliks <- matrix(NA, nrow(logliks), length(unique(merges)),
                       dimnames = list(rownames(logliks), unique(merges)))
  newlogliks <- sapply(unique(merges), function(newname) {
    oldnames <- names(merges)[merges == newname]
    newlogliks[, newname] = apply(logliks[, oldnames, drop = FALSE], 1, max)
  })
  newlogliks <- cbind(newlogliks, logliks[, setdiff(colnames(logliks), names(merges)), drop = FALSE])

  ## convert to probs:
  probs <- logliks2probs(newlogliks)
  clust <- colnames(probs)[apply(probs, 1, which.max)]
  names(clust) <- rownames(probs)
  out <- list(clust = clust, probs = round(probs, 2))  # (rounding probs to save memory)
  return(out)
}


# get a probabilities matrix from a logliks matrix
#' @export
probs2logliks <- function(probs) {
  return(log(probs))
}


