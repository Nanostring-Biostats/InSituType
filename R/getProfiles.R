
#' Extract mean background-subtracted profiles of RNA data
#'
#' Given cell assignments and count data, estimate the mean
#'  profile of each cluster.
#' 
#' @param x Counts matrix, cells * genes.
#' @param clust Vector of cluster assignments, or a matrix of probabilities
#'   of cells (rows) belonging to clusters (columns).
#' @param neg Vector of mean background counts (or a single value applied to all cells)
#' @return A matrix of gene x cell type expression profiles. 
#' @export
getRNAprofiles <- function(x, neg, clust) {
  if (length(neg) == 1) {
    neg <- rep(neg, nrow(x))
  }
  temp <- Estep(counts = x, clust = clust, neg = neg, assay_type = "RNA")
  return(temp$profiles)
}

#' Extract mean background-subtracted profiles of RNA data
#'
#' Given cell assignments and count data, estimate the mean
#'  profile of each cluster.
#' @param x Expression matrix, cells * proteins.
#' @param clust Vector of cluster assignments, or a matrix of probabilities
#'   of cells (rows) belonging to clusters (columns).
#' @param neg Vector of mean background counts
#' @return List with two elements: "profiles", a matrix of protein x cell type expression profiles, and "sds", a matrix of SD's.
#' @export
getProteinParameters <- function(x, clust) {
  temp <- Estep(counts = x, clust = clust, assay_type = "protein")
  return(temp)
}
