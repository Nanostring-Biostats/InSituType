% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/insitutypeML.R
\name{insitutypeML}
\alias{insitutypeML}
\alias{insitutypeML,ANY-method}
\alias{insitutypeML,SingleCellExperiment-method}
\title{Classify cells based on reference profiles}
\usage{
insitutypeML(x, ...)

\S4method{insitutypeML}{ANY}(
  x,
  neg = NULL,
  bg = NULL,
  cohort = NULL,
  reference_profiles,
  nb_size = 10,
  align_genes = TRUE
)

\S4method{insitutypeML}{SingleCellExperiment}(x, ..., assay.type = "counts")
}
\arguments{
\item{x}{Counts matrix (or dgCMatrix), cells * genes.

  Alternatively, a \linkS4class{SingleCellExperiment} object containing such
  a matrix.}

\item{...}{For the \linkS4class{SingleCellExperiment} method, additional
arguments to pass to the ANY method.}

\item{neg}{Vector of mean negprobe counts per cell. Can be provided}

\item{bg}{Expected background}

\item{cohort}{Vector of cells' cohort memberships}

\item{reference_profiles}{Matrix of expression profiles of pre-defined clusters,
e.g. from previous scRNA-seq. These profiles will not be updated by the EM algorithm.
Colnames must all be included in the init_clust variable.}

\item{nb_size}{The size parameter to assume for the NB distribution.}

\item{align_genes}{Logical, for whether to align the counts matrix and the reference_profiles by gene ID.}

\item{assay.type}{A string specifying which assay values to use.}
}
\value{
A list, with the following elements:
\enumerate{
\item clust: a vector given cells' cluster assignments
\item prob: a vector giving the confidence in each cell's cluster
\item profiles: Matrix of clusters' mean background-subtracted profiles
\item logliks: Matrix of cells' log-likelihoods under each cluster. Cells in rows, clusters in columns.
}
}
\description{
Supervised classification of cells. Each cell is assigned to the cell type 
 under which its observed expression profile is most likely.
}
\examples{
data("mini_nsclc")
data("ioprofiles")
sup <- insitutypeML(
 x = mini_nsclc$counts,
 neg = Matrix::rowMeans(mini_nsclc$neg),
 reference_profiles = ioprofiles)
table(sup$clust)
}
