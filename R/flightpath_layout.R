

#' "Flightpath" (umap-like) plot of clustering results
#'
#' Arrays cells in 2d space based on their probability of belonging to a given
#' cluster.
#' @param logliks Matrix of cells' log-likelihoods under each cluster. Must
#'   provide this or probs argument.
#' @param probs Matrix of cells' probabilities of belonging to each cluster.
#'   Must provide this or logliks argument.
#' @param profiles Matrix of cell type mean expression profiles. If provided,
#'   profiles rather than probs will be used to lay out the centroids.
#' @param cluster_xpos Vector of cluster centroids' x positions (i.e. where you
#'   want each cell type to appear in the plot)
#' @param cluster_ypos Vector of cluster centroids' y positions
#' @return A list with two elements: \enumerate{ \item clustpos: a matrix of
#'   cluster centroids * x,y positions in the flightpath plot \item cellpos: A
#'   matrix of cells * x,y positions in the flightpath plot }
#' @importFrom umap umap
#' @importFrom stats rnorm
#' @export
#' @examples
#' data("mini_nsclc")
#' unsup <- insitutype(
#'  x = mini_nsclc$counts,
#'  neg = Matrix::rowMeans(mini_nsclc$neg),
#'  assay_type = "RNA",
#'  n_clusts = 8,
#'  n_phase1 = 200,
#'  n_phase2 = 500,
#'  n_phase3 = 2000,
#'  n_starts = 1,
#'  max_iters = 5
#' ) # choosing inadvisably low numbers to speed the vignette; using the defaults in recommended.
#' flightpath_layout(logliks = unsup$logliks, profiles = unsup$profiles)
flightpath_layout <- function(logliks = NULL, probs = NULL, profiles = NULL, cluster_xpos = NULL, cluster_ypos = NULL) {

  if (is.null(probs) && is.null(logliks)) {
    stop("Must provide either probs or logliks.")
  }
  if (is.null(probs) && !is.null(logliks)) {
    probs <- logliks2probs(logliks)
  }
  # force NA probs to 0:
  probs <- replace(probs, is.na(probs), 0)
  # get cluster centroid positions if not pre-specified:
  if (is.null(cluster_xpos) || is.null(cluster_ypos)) {
    # controls for a umap-based layout:
    conf <- umap::umap.defaults
    conf$min_dist <- 3
    conf$spread <- conf$min_dist * 1.1
    conf$n_neighbors <- ncol(probs)
    if (!is.null(profiles)) {
      clustum <- umap::umap(t(sqrt(profiles)), config = conf)$layout
    } else {
      clustum <- umap::umap(t(probs), config = conf)$layout
    }
    
    cluster_xpos <- clustum[, 1]
    cluster_ypos <- clustum[, 2]
  }

  # get cell xy positions as a weighted average of the umap positions
  ux <- probs %*% cluster_xpos
  uy <- probs %*% cluster_ypos
  
  # jitter the xy positions, jittering widely for prob = 1 cells and minimally for prob < 0.5 cells:
  jitterrange <- 0.01 * c(0.0005, 0.9) * max(diff(range(ux)), diff(range(uy))) 
  jitteramount <- jitterrange[1] + pmax((2 * apply(probs, 1, max) - 1), 0)  * jitterrange[2]
  ux <- ux + rnorm(length(ux), mean = 0, sd = jitteramount)
  uy <- uy + rnorm(length(ux), mean = 0, sd = jitteramount)

  out <- list(clustpos = cbind(cluster_xpos, cluster_ypos),
             cellpos = cbind(ux, uy),
             clust = colnames(probs)[apply(probs, 1, which.max)])
  colnames(out$clustpos) <- c("x", "y")
  colnames(out$cellpos) <- c("x", "y")
  
  # get clusters' mean confidence:
  out$meanconfidence <- getMeanClusterConfidence(probs)
  return(out)
}




#'Plot flightpath results
#'
#'@param flightpath_result The list output by the flightpath_layout function.
#'  Two elements: clustpos, cellpos. Must provide either this or
#'  insitutype_result.
#'@param insitutype_result The list output by insitutype or insitutypeML. Must
#'  provide either this or insitutype_result.
#'@param col Optional, a vector of cell colors, with length equal to the number
#'  of individual cells.
#'@param showclusterconfidence Logical, for whether to label clusters with the
#'  average posterior probability of the cells within them. Gives a readout of
#'  how distinct a cluster is from the others.
#'@importFrom utils data
#'@importFrom scales alpha
#'@import ggplot2
#'@importFrom grDevices colors
#'@importFrom rlang .data
#'@return a ggplot object
#'
#'@export
#'@examples 
#' data("ioprofiles")
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
#' flightpath_plot(insitutype_result = unsup)

flightpath_plot <- function(flightpath_result = NULL, insitutype_result = NULL, col = NULL, showclusterconfidence = TRUE){
  
  # get the flightpath results to use 
  if (!is.null(flightpath_result) && !is.null(insitutype_result)) {
    warning("flightpath_result and insitutype_result were both provided. Using only flightpath_result.")
    insitutype_result <- NULL
  }
  if (is.null(flightpath_result) && is.null(insitutype_result)) {
    stop("Must provide either flightpath_result or insitutype_result.")
  }
  if (is.null(flightpath_result)) {
    flightpath_result <- flightpath_layout(logliks = insitutype_result$logliks, profiles = insitutype_result$profiles)
  }
  
  # create color scheme if needed:
  if (is.null(col)) {
    utils::data("iocolors", package = "InSituType", envir = environment())
    scols <- c('#8DD3C7','#FFFFB3','#BEBADA','#FB8072','#80B1D3','#FDB462','#B3DE69','#FCCDE5','#D9D9D9','#BC80BD',
               '#CCEBC5','#FFED6F','#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF','#999999', 
               sample(colors()[!grepl("grey", colors())], 100))[seq_along(unique(flightpath_result$clust))]
    names(scols) <- unique(flightpath_result$clust)
    iotypespresent <- intersect(names(environment()[['iocolors']]), names(scols))
    scols[iotypespresent] <- environment()[['iocolors']][iotypespresent]
    col <- scols[flightpath_result$clust]
  }

  # prep data for plotting:
  df <-
    data.frame(
      x = flightpath_result$cellpos[, 1],
      y = flightpath_result$cellpos[, 2],
      col = scales::alpha(col, 0.7)
    )
  df_text <- data.frame(x = flightpath_result$clustpos[, 1],
                        y = flightpath_result$clustpos[, 2],
                        group = rownames(flightpath_result$clustpos),
                        col = "black")
  
  if (showclusterconfidence) {
    confthresh <- 0.8
    confidencecolors <- c('#FEB24C','#FD9D43','#FC863A','#FC6330','#F64226',
                          '#E8251F','#D2111F','#B60224','#620015','#000000')
    df_text$col <- confidencecolors[
      1 + round(9 * (pmax(flightpath_result$meanconfidence, confthresh) - confthresh) / (1 - confthresh))]
    
    df_text$group <- paste0(df_text$group, "(", round(flightpath_result$meanconfidence, 2), ")")
  }
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(df, mapping  = ggplot2::aes(x = flightpath_result$cellpos[, 1], 
                                                    y = flightpath_result$cellpos[, 2], 
                                                    color = I(col),
                                                    size = I(0.1))) +
    ggplot2::scale_color_identity() +
    ggplot2::geom_text(df_text,
              mapping = ggplot2::aes(x = .data$x, y = .data$y, label = .data$group, col = I(col)),
              size = 3) +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
          panel.grid = ggplot2::element_blank(),
          axis.text = ggplot2::element_blank())
  flightpath_plot_folder <- "./NBClust-Plots" # tempdir()
  if (!dir.exists(flightpath_plot_folder)) dir.create(flightpath_plot_folder, showWarnings = FALSE, recursive = TRUE)
  flightpath_plot_filename <- paste(format(Sys.time(), "%Y-%m-%d_%H-%M-%S-%Z"), "flightpath_plot.png", sep="-")
  flightpath_plot_file <- paste(flightpath_plot_folder,flightpath_plot_filename , sep="/")
  message("Saving flightpath_plot to: ", flightpath_plot_file)
  ggsave(filename = flightpath_plot_filename, plot = p, device = "png", path = flightpath_plot_folder,
         width = 7,
         height = 7,
         units="in")

  return(p)
}


#' Summarize clusters' mean confidence
#' 
#' Calculate the mean confidence of the cell calls from each cluster
#' @param probs Matrix of probabilities
#' @return a vector of mean confidences, with values of 1 corresponding to clusters with only prob == 1
#' @examples
#' data("mini_nsclc")
#' probs <- sapply(rownames(mini_nsclc$counts), function(x) {a = runif(10); a/sum(a)})
#' dimnames(probs)[[1]] <- letters[1:10]
#' probs <- t(probs)
#' getMeanClusterConfidence(probs)
getMeanClusterConfidence <- function(probs) {
  
  maxprobs <- apply(probs, 1, max, na.rm = TRUE)
  meanconfidence <- sapply(colnames(probs), function(name) {
    thisclust <- probs[, name] == maxprobs
    mean(probs[thisclust, name, drop = FALSE])
  })
  
  return(meanconfidence)
}
