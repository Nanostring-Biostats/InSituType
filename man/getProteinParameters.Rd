% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getProfiles.R
\name{getProteinParameters}
\alias{getProteinParameters}
\title{Extract mean background-subtracted profiles of RNA data}
\usage{
getProteinParameters(x, clust)
}
\arguments{
\item{x}{Expression matrix, cells * proteins.}

\item{clust}{Vector of cluster assignments, or a matrix of probabilities
of cells (rows) belonging to clusters (columns).}

\item{neg}{Vector of mean background counts}
}
\value{
List with two elements: "profiles", a matrix of protein x cell type expression profiles, and "sds", a matrix of SD's.
}
\description{
Given cell assignments and count data, estimate the mean
 profile of each cluster.
}
