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
#' @return A profiles matrix with the rows rescaled according to platform effects, 
#'  and a vector of genes' scaling factors (what they were multiplied by when updating the reference profiles). 
#' @export
#' @importFrom MASS glm.nb
rescaleProfiles <- function(counts, neg, fixed_profiles, align_genes = TRUE, nb_size = 10) {
  
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
  res <- spatialdecon(norm = cbind(totcounts, totcounts), 
                      bg = sum(neg), 
                      X = fixed_profiles, 
                      resid_thresh = Inf)
  scaling_factors <- 2^res$resids[, 1]
  rescaledprofiles <- sweep(fixed_profiles, 1, scaling_factors, "*")
  
  out = list(profiles = rescaledprofiles,
             scaling_factors = scaling_factors)
  return(out)
}
  