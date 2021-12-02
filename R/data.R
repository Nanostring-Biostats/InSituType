
#' Small example SMI data from a NSCLC tumor
#'
#' A 2000-cell excerpt from a 1000-plex SMI study of a NSCLC tumor. 
#'
#' @format A list with the following elements:
#'  \itemize{
#'  \item counts A matrix of raw counts, with cells in rows and genes in columns
#'  \item counts A matrix of negprobe counts, with cells in rows and negprobes in columns
#'  \item x x positions
#'  \item y y position
#'  \item umap umap projection
#'  }
"mini_nsclc"



#' Matrix of immune cell profiles
#'
#' A matrix of gene * cell type expected expression values
#'
#' @format A matrix of 27161 genes x 16 cell types. 
"ioprofiles"

#' Default colors for the cell types in the ioprofiles matrix
#'
#' A named vector of colors, giving colors for the cell types of the ioprofiles
#'  matrix.
#'
#' @format A named vector
"iocolors"