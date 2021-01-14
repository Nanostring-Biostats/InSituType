

#' "Flightpath" (umap-like) plot of clustering results
#' 
#' Arrays cells in 2d space based on their probability of belonging to a given cluster. 
#' @param probs Matrix of cells' probabilities of belonging to each cluster.
#' @param cluster_xpos Vector of cluster centroids' x position (i.e. where you want each cell type to appear in the plot)
#' @param cluster_ypos Vector of cluster centroids' y position
#' @return A list with two elements:
#' \enumerate{
#' \item clustpos: a matrix of cluster centroids * x,y positions in the flightpath plot
#' \item cellpos: A matrix of cells * x,y positions in the flightpath plot
#' }
#' @importFrom umap umap
#' @export
flightpath_layout <- function(probs, cluster_xpos = NULL, cluster_ypos = NULL) {
  
  # get cluster centroid positions if not pre-specified:
  if (is.null(cluster_xpos) | is.null(cluster_ypos)) {
    # controls for a umap-based layout:
    conf = umap.defaults
    conf$min_dist = 3
    conf$spread = conf$min_dist * 1.1
    conf$n_neighbors = ncol(semi$profiles)
    clustum = umap(t(sqrt(semi$profiles)), config = conf)$layout
    cluster_xpos = clustum[, 1]
    cluster_ypos = clustum[, 2]
  }
  
  # get cell xy positions as a weighted average of the umap positions
  ux = semi$probs %*% clustum[, 1]
  uy = semi$probs %*% clustum[, 2]
  
  out = list(clustpos = cbind(cluster_xpos, cluster_ypos),
             cellpos = cbind(ux, uy))
  colnames(out$clustpos) = c("x", "y")
  colnames(out$cellpos) = c("x", "y")
  return(out)
}



#'Generate a probability-based umap from the cell typing results 
#'
#'@param semi the cell typing result from the cellEMClust function
#'
#'@importFrom umap umap.defaults
#'@importFrom umap umap
#'@importFrom SpatialDecon cellcols
#'@importFrom scales alpha
#'@importFrom RColorBrewer brewer.pal
#'
#'@return a ggplot object
#'
#'@export
#'
plotCellTypeUmap <- function(semi){
  # create color schemes in addition to cellcols
  cellcols <- SpatialDecon::cellcols
  scols <- c(cellcols, RColorBrewer::brewer.pal(12, "Set3")[1:length(setdiff(unique(semi$clust), names(cellcols)))])
  names(scols) <- c(names(cellcols), setdiff(unique(semi$clust), names(cellcols)))
  
  # change umap default settings 
  par(mfrow <- c(1, 1))
  conf <- umap::umap.defaults
  conf$min_dist <- 3
  conf$spread = conf$min_dist * 1.1
  conf$n_neighbors <- ncol(semi$profiles)
  #plot(umap(t(semi$profiles), config = conf)$layout)
  clustum <- umap::umap(t(sqrt(semi$profiles)), config = conf)$layout
  
  # get cell xy positions as a weighted average of the umap positions
  ux <- semi$probs %*% clustum[, 1]
  uy <- semi$probs %*% clustum[, 2]
  
  df <- data.frame(ux, uy, col = scales::alpha(scols[semi$clust], 0.7))
  df_text <- data.frame(x = clustum[, 1], 
                        y = clustum[, 2], 
                        group = colnames(semi$probs))
  p <- ggplot() +
    geom_point(df, mapping  = aes(x = ux, y = uy, color = col)) +
    scale_color_identity() +
    geom_text(df_text, 
              mapping = aes(x = x, y = y, label = group),
              size = 5) +
    xlab("") +
    ylab("") +
    theme_bw() +
    theme(legend.position = "none", 
          panel.grid = element_blank(),
          axis.text = element_blank())
  
  return(p)
}