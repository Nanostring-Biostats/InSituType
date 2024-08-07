#' Get the neighborhood expression profile around all cells
#' @param counts Counts matrix
#' @param xy 2-column matrix of cells' xy positions
#' @param tissue vector of tissue IDs. Used to ensure cells for different tissues are never called neighbors
#' @param N number of neighbors to use. Specify this or \code{rad}. 
#' @param rad radius to use to define neighbors. Specify this or \code{N}. 
#' @param dim_reduce_to If entered, the neighborhood matrix will be reduced to this many PCs
#' @return A matrix of neighborhood expression, potentially by gene, or else by PCs if \code{dim_reduce_to} was set.
#' @export
#' @importFrom irlba prcomp_irlba
getSpatialContext <- function(counts, xy, tissue = NULL, N = 50, rad = NULL, dim_reduce_to = NULL) {
  
  # define neighbors:
  if (is.null(tissue)) {
    tissue = 1
  }
  if (!is.null(N)) {
    neighbors <- nearestNeighborGraph(x = xy[, 1], y = xy[, 2], N = N, subset = tissue) 
    rad <- NULL
  } 
  if (!is.null(rad)) {
    neighbors <- radiusBasedGraph(x = xy[, 1], y = xy[, 2], R = rad, subset = tissue) 
  } 
  
  # get neighborhood expression:
  neighborexpression <- get_neighborhood_expression(counts = counts, neighbors = neighbors) 
  
  # dimension reduce
  if (!is.null(dim_reduce_to)) {
    neighborexpression <- irlba::prcomp_irlba(neighborexpression, n = dim_reduce_to)$x
  }
  return(neighborexpression)
}




#' Create spatial network from N nearest neighbors
#'
#' For each cell identify \code{N} nearest neighbors in Euclidean space and
#' create an edge between them in graph structure, optionally subset cells (see
#' Details).
#'
#' Edges will only be created for cells that have the same \code{subset} value,
#' usually the slide column id but could also be a slide plus FOV id to only
#' create edges within an FOV.
#'
#' @param x spatial coordinate
#' @param y spatial coordinate
#' @param N number of nearest neighbors
#' @param subset same length as x,y (see Details)
#'
#' @return sparse adjacency matrix with distances
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
#' @importFrom spatstat.geom nnwhich
#' @importFrom spatstat.geom nndist
#' @importFrom Matrix sparseMatrix
nearestNeighborGraph <- function(x, y, N, subset=1) {
  DT <- data.table::data.table(x = x, y = y, subset = subset)
  nearestNeighbor <- function(i) {
    subset_dt <- DT[subset == i]
    idx <- which(DT[["subset"]] == i)
    ndist <- spatstat.geom::nndist(subset_dt[, .(x, y)],
                                   k=1:N)
    nwhich <- spatstat.geom::nnwhich(subset_dt[, .(x, y)],
                                     k=1:N)
    ij <- data.table::data.table(i = idx[1:nrow(subset_dt)],
                                 j = idx[as.vector(nwhich)],
                                 x = as.vector(ndist))
    return(ij)
  }
  ij <- data.table::rbindlist(lapply(unique(subset), nearestNeighbor))
  adj.m <- Matrix::sparseMatrix(i = ij$i, j = ij$j, x = ij$x, dims = c(nrow(DT), nrow(DT)))
  return(adj.m)
}

#' Create spatial network from neighbors within radius R
#'
#' For each cell identify neighbors within distance \code{R} in Euclidean space
#' and create an edge between them in graph structure, optionally subset cells
#' (see Details).
#'
#' Edges will only be created for cells that have the same \code{subset} value,
#' usually the slide column id but could also be a slide plus FOV id to only
#' create edges within an FOV.
#'
#' @param x spatial coordinate
#' @param y spatial coordinate
#' @param R radius
#' @param subset same length as x,y (see Details)
#'
#' @return sparse adjacency matrix with distances
#' @importFrom data.table data.table
#' @importFrom data.table rbindlist
#' @importFrom Matrix sparseMatrix
#' @importFrom spatstat.geom ppp
#' @importFrom spatstat.geom closepairs
radiusBasedGraph <- function(x, y, R, subset=1) {
  DT <- data.table::data.table(x = x, y = y, subset = subset)
  radiusNeighbor <- function(i) {
    subset_dt <- DT[subset == i]
    idx <- which(DT[["subset"]] == i)
    pp <- spatstat.geom::ppp(subset_dt$x, subset_dt$y,
                             range(subset_dt$x), range(subset_dt$y))
    cp <- spatstat.geom::closepairs(pp, R)
    ij <- data.table::data.table(i = idx[cp$i],
                                 j = idx[cp$j],
                                 x = cp$d)
    return(ij)
  }
  ij <- data.table::rbindlist(lapply(unique(subset), radiusNeighbor))
  adj.m <- Matrix::sparseMatrix(i = ij$i, j = ij$j, x = ij$x, dims = c(nrow(DT), nrow(DT)))
  return(adj.m)
}



#' Calculate neighborhood expression
#'
#' Calculates the expression profile of each cell's neighborhood
#' @param counts Single cell expression matrix
#' @param neighbors A neighbors adjacency matrix
#' @return A matrix in the same dimensions as \code{counts}, giving the expression profile of each cell's neighborhood.
get_neighborhood_expression <- function(counts, neighbors) {
  
  # check:
  if (nrow(counts) != ncol(neighbors)) {
    stop("misalignment between nrow(counts) and ncol(neighbors)")
  }
  # get clust-specific environment expression
  env <- neighbor_colMeans(counts, neighbors)
  rownames(env) <- rownames(neighbors)
  env <- as.matrix(env)
  return(env)
}

#' for each cell, get the colMeans of x over its neighbors:
#' @param x A matrix
#' @param neighbors A (probably sparse) adjacency matrix
neighbor_colMeans <- function(x, neighbors) {
  neighbors@x <- rep(1, length(neighbors@x))
  neighbors <- Matrix::Diagonal(x=1/Matrix::rowSums(neighbors)) %*% neighbors
  neighbors@x[neighbors@x==0] <- 1
  out <- neighbors %*% x
  return(out)
}



#' for each cell, get the colSums of x over its neighbors:
#' @param x A matrix
#' @param neighbors A (probably sparse) adjacency matrix
neighbor_colSums <- function(x, neighbors) {
  neighbors@x <- rep(1, length(neighbors@x))
  neighbors <- Matrix::Diagonal(x=rep(1, nrow(neighbors))) %*% neighbors
  neighbors@x[neighbors@x==0] <- 1
  out <- neighbors %*% x
  return(out)
}
