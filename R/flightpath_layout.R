

#' "Flightpath" (umap-like) plot of clustering results
#'
#' Arrays cells in 2d space based on their probability of belonging to a given cluster.
#' @param probs Matrix of cells' probabilities of belonging to each cluster.
#' @param profiles Matrix of cell type mean expression profiles. If provided, profiles rather than probs will be used to lay out the centroids. 
#' @param cluster_xpos Vector of cluster centroids' x positions (i.e. where you want each cell type to appear in the plot)
#' @param cluster_ypos Vector of cluster centroids' y positions
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
    # controls for a umap-based layout:
    conf = umap.defaults
    conf$min_dist = 3
    conf$spread = conf$min_dist * 1.1
    conf$n_neighbors = ncol(probs)
    if (!is.null(profiles)) {
      clustum = umap(t(sqrt(profiles)), config = conf)$layout
    } else {
      clustum = umap(t(probs), config = conf)$layout
    }
    
    cluster_xpos = clustum[, 1]
    cluster_ypos = clustum[, 2]
  }

  # get cell xy positions as a weighted average of the umap positions
  ux = probs %*% clustum[, 1]
  uy = probs %*% clustum[, 2]

  out = list(clustpos = cbind(cluster_xpos, cluster_ypos),
             cellpos = cbind(ux, uy),
             clust = colnames(probs)[apply(probs, 1, which.max)])
  colnames(out$clustpos) = c("x", "y")
  colnames(out$cellpos) = c("x", "y")
  return(out)
}




#' Plot flightpath results
#'
#'@param fp The list output by the flightpath_layout function. Two elements: clustpos, cellpos and clust
#'@param col Optional, a vector of cell colors 
#'@importFrom utils data
#'@importFrom scales alpha
#'@importFrom RColorBrewer brewer.pal
#'
#'@return a ggplot object
#'
#'@export
#'
flightpath_plot <- function(flightpath_result = NULL, insitutype_result = NULL, col = NULL){
  
  # get the flightpath results to use 
  if (!is.null(flightpath_result) & !is.null(insitutype_result)) {
    warning("flightpath_result and insitutype_result were both provided. Using only flightpath_result.")
    insitutype_result <- NULL
  }
  if (is.null(flightpath_result) & is.null(insitutype_result)) {
    stop("Must provide either flightpath_result or insitutype_result.")
  }
  if (is.null(flightpath_result)) {
    flightpath_result <- flightpath_layout(probs = res$probs, profiles = res$profiles)
  }
  
  # create color scheme if needed:
  if (is.null(col)) {
    utils::data("iocolors", envir = environment())
    scols <- c(RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(9, "Set1"), sample(colors()[!grepl("grey", colors())], 100))[1:length(unique(flightpath_result$clust))]
    names(scols) <- unique(flightpath_result$clust)
    iotypespresent = intersect(names(iocolors), names(scols))
    scols[iotypespresent] = iocolors[iotypespresent]
    col = scols[flightpath_result$clust]
  }

  # prep data for plotting:
  df <- data.frame(x = flightpath_result$cellpos[, 1], y = flightpath_result$cellpos[, 2], col = scales::alpha(col, 0.7))
  df_text <- data.frame(x = flightpath_result$clustpos[, 1],
                        y = flightpath_result$clustpos[, 2],
                        group = rownames(flightpath_result$clustpos))
  p <- ggplot() +
    geom_point(df, mapping  = aes(x = flightpath_result$cellpos[, 1], y = flightpath_result$cellpos[, 2], color = I(col))) +
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


