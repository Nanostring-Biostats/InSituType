#' Rescale rows of reference profiles
#' 
#' Estimates platform effects and rescales the reference matrix accordingly
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell. Can be provided 
#' @param fixed_profiles Matrix of expression profiles of pre-defined clusters,
#'  e.g. from previous scRNA-seq. These profiles will not be updated by the EM algorithm.
#'  Colnames must all be included in the init_clust variable.
#' @param align_genes Logical, for whether to align the counts matrix and the fixed_profiles by gene ID.
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param max_rescaling Scaling factors will be truncated above by this value and below by its inverse (at 1/value and value)
#' @return 
#' \enumerate{
#' \item profiles: A profiles matrix with the rows rescaled according to platform effects
#' \item scaling_factors: A vector of genes' scaling factors (what they were multiplied by when updating the reference profiles). 
#' }
#' @export
#' @importFrom SpatialDecon spatialdecon
rescaleProfiles <- function(counts, neg, fixed_profiles, align_genes = TRUE, nb_size = 10, max_rescaling = 5) {
  
  # align genes:
  if (align_genes) {
    sharedgenes <- intersect(rownames(fixed_profiles), colnames(counts))
    lostgenes <- setdiff(colnames(counts), rownames(fixed_profiles))
    
    # subset:
    counts <- counts[, sharedgenes]
    fixed_profiles <- fixed_profiles[sharedgenes, ]
    
    # warn about genes being lost:
    if ((length(lostgenes) > 0) & length(lostgenes < 50)) {
      message(paste0("The following genes in the count data are missing from fixed_profiles and will be omitted from anchor selection: ",
                     paste0(lostgenes, collapse = ",")))
    }
    if (length(lostgenes) > 50) {
      message(paste0(length(lostgenes), " genes in the count data are missing from fixed_profiles and will be omitted from anchor selection"))
    }
  }
  
  # get background-subtracted bulk profile:
  totcounts <- matrix(colSums(counts), ncol(counts))
  rownames(totcounts) <- colnames(counts)
  totcountsbgsub <- pmax(totcounts - sum(neg), 0)
  
  # deconvolve the bulk profile with spatialdecon:
  res <- SpatialDecon::spatialdecon(norm = cbind(totcounts, totcounts), 
                      bg = sum(neg), 
                      X = fixed_profiles, 
                      resid_thresh = Inf)
  log2resids <- res$resids[, 1]
  log2resids <- log2resids - mean(log2resids)
  scaling_factors <- 2^log2resids
  scaling_factors <- pmax(pmin(scaling_factors, max_rescaling), max_rescaling^-1)
  rescaledprofiles <- sweep(fixed_profiles, 1, scaling_factors, "*")
  
  out = list(profiles = rescaledprofiles,
             scaling_factors = scaling_factors)
  return(out)
}
  



#' Use anchor cells to update reference profiles
#' 
#' Uses anchor cells to estimate platform effects / scaling factors to be applied to
#'  the genes/rows of the reference profile matrix. Then uses Bayesian math to
#'  update the individual elements on X. 
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell. Can be provided 
#' @param anchors Vector of anchor assignments
#' @param fixed_profiles Matrix of expression profiles of pre-defined clusters,
#'  e.g. from previous scRNA-seq. These profiles will not be updated by the EM algorithm.
#'  Colnames must all be included in the init_clust variable.
#' @param align_genes Logical, for whether to align the counts matrix and the fixed_profiles by gene ID.
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param max_rescaling Scaling factors will be truncated above by this value and below by its inverse (at 1/value and value)
#' @return 
#' \enumerate{
#' \item profiles: A profiles matrix with the rows rescaled according to platform effects and individual elements updated further
#' \item scaling_factors: A vector of genes' scaling factors (what they were multiplied by when updating the reference profiles). 
#' }
#' @export
rescaleProfiles <- function(counts, neg, anchors, fixed_profiles, align_genes = TRUE, nb_size = 10, max_rescaling = 5) {
  
  # get total expression and negs per anchor cell type:
  sums <- sapply(by(counts, anchors, colSums), cbind)
  rownames(sums) <- colnames(counts)
  temp <- by(neg, anchors, sum)
  negsums <- as.vector(temp)
  names(negsums) <- names(temp)
  negsums <- negsums[colnames(sums)]
  
  ### build data frame for estimating row/gene scaling factors:
  celltypes <- intersect(colnames(fixed_profiles), unique(anchors))
  genes <- intersect(rownames(sums), rownames(fixed_profiles))
  df <- data.frame(
    gene = rep(genes, length(celltypes)),
    celltype = rep(celltypes, each = lenght(genes)),
    ref = as.vector(fixed_profiles[genes, celltypes]),
    rnasum = as.vector(sums[genes, celltypes]),
    negsum = rep(negsums, each = length(genes))
  )
  df$bgsub <- pmax(df$rnasum - df$negsum, 0)
  # throw out rows with values too low to be useful in estimating a log-ratio:
  
  
  
  
  
}