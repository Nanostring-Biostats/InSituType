#' "Flightpath" (umap-like) plot of clustering results
#' 
#' Arrays cells in 2d space based on their probability of belonging to a given cluster. 
#' @param probs Matrix of cells' probabilities of belonging to each cluster.
#' @param profiles Matrix of mean cluster profiles, used in umap layout. 
#'   If not provided, then the probs matrix is passed to umap.
#' @param cluster_xpos Vector of cluster centroids' x position (i.e. where you want each cell type to appear in the plot)
#' @param cluster_ypos Vector of cluster centroids' y position
#' @return A list with two elements:
#' \enumerate{
#' \item clustpos: a matrix of cluster centroids * x,y positions in the flightpath plot
#' \item cellpos: A matrix of cells * x,y positions in the flightpath plot
#' }
#' @importFrom umap umap
#' @export
flightpath_layout <- function(probs, profiles = NULL, cluster_xpos = NULL, cluster_ypos = NULL) {
  
  
  # get cluster centroid positions if not pre-specified:
  if (is.null(cluster_xpos) | is.null(cluster_ypos)) {
    # matrix for umap layout:
    mat = probs
    if (!is.null(profiles)) {
      mat = profiles
    }
    # controls for a umap-based layout:
    conf = umap::umap.defaults
    conf$min_dist = 3
    conf$spread = conf$min_dist * 1.1
    conf$n_neighbors = ncol(mat)
    clustum = umap(t(sqrt(mat)), config = conf)$layout
    cluster_xpos = clustum[, 1]
    cluster_ypos = clustum[, 2]
  }
  
  # get cell xy positions as a weighted average of the umap positions
  ux = probs %*% clustum[, 1]
  uy = probs %*% clustum[, 2]
  
  out = list(clustpos = cbind(cluster_xpos, cluster_ypos),
             cellpos = cbind(ux, uy))
  colnames(out$clustpos) = c("x", "y")
  colnames(out$cellpos) = c("x", "y")
  return(out)
}
