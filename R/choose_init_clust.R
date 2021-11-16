
#' Choose init clusts aware of fixed_profiles
#'
#' Choose init clusts such that cells with good logliks under fixed_profiles aren't initialized to new clusters 
#' @param counts Counts matrix from inside insitutype
#' @param fixed_profiles Matrix of expression profiles of pre-defined clusters,
#'  e.g. from previous scRNA-seq. These profiles will not be updated by the EM algorithm.
#'  Colnames must all be included in the init_clust variable.
#' @param bg Expected background
#' @param size The size parameter to assume for the NB distribution.
#' @param n_clusts Number of free clusters
#' @param thresh Threshold above which cells will be called as well-fit by fixed_profiles and initialized as a fixed_profile.
#' @return A vector of initial cluster assignments
choose_init_clust <- function(counts, fixed_profiles, bg, size, n_clusts, thresh = 0.9) {#<------------------------- consider thresh carefully
  
  if (align_genes & !is.null(fixed_profiles)) {
    sharedgenes <- intersect(rownames(fixed_profiles), colnames(counts))
    lostgenes <- setdiff(colnames(counts), rownames(fixed_profiles))
    
    # subset:
    counts <- counts[, sharedgenes]
    fixed_profiles <- fixed_profiles[sharedgenes, ]
    
  }
  # get logliks under all fixed profiles
  logliks <- Mstep(counts = counts, 
                   means = fixed_profiles, 
                   bg = bg, 
                   size = size, 
                   digits = 2, 
                   return_loglik = TRUE) 
    
  # determine a cutoff based on max loglik, tot counts
  maxloglik <- apply(logliks, 1, max)
  totcounts <- rowSums(counts)
  mod <- lm(maxloglik ~ totcounts)$coef
  #stat <- maxloglik / totcounts   #<------------------------- consider this carefully
  #assigntofixed <- stat > thresh
  assigntofixed <- maxloglik > mod[1] + thresh * mod[2] * totcounts
    
  # ensure that at least a reasonable number of cells are assigned to free clusters:
  
  # assign well-fit cells to fixed profiles
  init_clust <- rep(NA, nrow(counts))
  names(init_clust) <- rownames(counts)
  init_clust[assigntofixed] <- colnames(logliks)[apply(logliks[assigntofixed, ], 1, which.max)]
  
  # randomly initialize poorly-fit cells to new clusters
  newandfixedclusternames <- makeClusterNames( colnames(fixed_profiles) , ncol(fixed_profiles) + n_clusts)
  newclusternames <- setdiff(newandfixedclusternames, colnames(fixed_profiles))
  init_clust[!assigntofixed] <- sample(newclusternames, sum(!assigntofixed), replace = TRUE)
  
  return(init_clust)  
}
  

