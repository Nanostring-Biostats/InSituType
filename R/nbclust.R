
#' Calculate the likelihood of the vector x using the reference vector of y
#'
#' @param x a vector of a reference cell type
#' @param mat a matrix of expression levels in all cells
#' @param bg background level (default: 0.01)
#' @param size the parameters for dnbinom function (default: 10)
#' @param digits the number of digits for rounding
#'
#' @importFrom Matrix rowSums
#' @importFrom stats dnbinom
#'
#' @export
lldist <- function(x, mat, bg = 0.01, size = 10, digits = 2) {
  # convert to matrix form if only a vector was input:
  if (is.vector(mat)) {
    mat <- matrix(mat, nrow = 1)
  }
  # Check dimensions on bg and stop with informative error if not
  #  conformant
  if ( is.vector( bg ) )
  {
    if ( !identical( length( bg ) , nrow( mat ) ) )
    {
      errorMessage <- sprintf( "Dimensions of count matrix and background are not conformant.\nCount matrix rows: %d, length of bg: %d" , nrow( mat ) , length( bg ) )
      stop( errorMessage )
    }
  }

  # calc scaling factor to put y on the scale of x:
  if ( is.vector( bg ) )
  {
    bgsub <- pmax( sweep( mat , 1 , bg , "-" ) , 0 )
  }
  else
  {
    bgsub <- pmax( mat - bg , 0 )
  }
  s <- rowSums( bgsub ) / sum( x )
  # override it if s is negative:
  s[s <= 0] = Matrix::rowSums(mat[s <= 0, , drop = FALSE]) / sum(x)

  # expected counts:
  if ( is.vector( bg ) )
  {
    yhat <- sweep( s %*% t( x ) , 1 , bg , "+" )
  }
  else
  {
    yhat <- s %*% t(x) + bg
  }

  # loglik:
  lls <- stats::dnbinom(x = as.matrix(mat), size = size, mu = yhat, log = TRUE)

  return(round(rowSums(lls), digits))
}


#' M step
#'
#' Compute probability that each cell belongs to a given cluster
#' @param counts Counts matrix, cells * genes.
#' @param means Matrix of mean cluster profiles,
#'  with genes in rows and clusters in columns.
#' @param cohort a vector of cells' "cohort" assignment, used to update logliks 
#'  based on cluster frequencies within a cohort.
#' @param bg Expected background
#' @param size NB size parameter
#' @param digits Round the output to this many digits (saves memory)
#' @param return_loglik If TRUE, logliks will be returned. If FALSE, probabilities will be returned. 
#' @return Matrix of probabilities of each cell belonging to each cluster
#' @export
Mstep <- function(counts, means, cohort, bg = 0.01, size = 10, digits = 2, return_loglik = FALSE) {
  # get logliks of cells * clusters
  logliks <- apply(means, 2, function(x) {
    lldist(x = x, mat = counts, bg = bg, size = size)
  })
  # adjust by cohort frequency:
  logliks <- update_logliks_with_cohort_freqs(logliks = logliks, 
                                              cohort = cohort, 
                                              minfreq = 1e-4, 
                                              nbaselinecells = 100) 
  if (return_loglik) {
    return(round(logliks, digits))
  } else {
    # first rescale (ie recenter on log scale) to avoid rounding errors:
    logliks <- sweep( logliks , 1 , apply( logliks , 1 , max ) , "-" )
    # get on likelihood scale:
    liks <- exp(logliks)
    # convert to probs
    probs <- sweep( liks , 1 , rowSums( liks ) , "/" )
    return(round(probs, digits))
  }
}


#' E step: estimate each cluster's mean profile
#'
#' Given cell assignments (or posterior probabilities), estimate the mean
#'  profile of each cluster.
#' @param counts Counts matrix, cells * genes.
#' @param clust Vector of cluster assignments, or a matrix of probabilities
#'   of cells (rows) belonging to clusters (columns).
#' @param neg Vector of mean background counts
#'
#' @importFrom Matrix rowSums
#'
#' @return A matrix of cluster profiles, genes * clusters
#' @export
Estep <- function(counts, clust, neg) {

  # get cluster means:
  means <- sapply(unique(clust), function(cl) {
    pmax(Matrix::colMeans(counts[clust == cl, , drop = FALSE]) - mean(neg[clust == cl]), 0)
  })
  return(means)
}





#' Cluster via EM algorithm based on cell logliks
#'
#' Cluster single cell gene expression data using an EM algorithm.
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param bg Expected background
#' @param cohort Vector of cells' "cohort" assignments, uses to assess frequencies in each cluster. 
#' @param init_profiles Matrix of cluster profiles under which to begin iterations.
#' If NULL, initial assignments will be automatically inferred, using init_clust 
#' if available, and using random clusters if not. 
#' @param init_clust Vector of initial cluster assignments.
#' If NULL, initial assignments will be automatically inferred.
#' @param n_clusts Number of clusters, in addition to any pre-specified cell types.
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param  Numer of iterations
#' @param pct_drop the decrease in percentage of cell types with a valid switchover to 
#'  another cell type compared to the last iteration. Default value: 1/10000. A valid 
#'  switchover is only applicable when a cell has changed the assigned cell type with its
#'  highest cell type probability increased by min_prob_increase. 
#' @param min_prob_increase the threshold of probability used to determine a valid cell 
#'  type switchover
#' @param max_iters Maximum number of iterations
#'
#'  @importFrom stats lm
#'
#' @return A list, with the following elements:
#' \enumerate{
#' \item probs: a matrix of probabilities of all cells (rows) belonging to all clusters (columns)
#' \item profiles: a matrix of cluster-specific expression profiles
#' }
#' @export
nbclust <- function(counts, neg, bg = NULL, 
                    fixed_profiles = NULL,
                    init_profiles = NULL, init_clust = NULL, n_clusts = NULL,
                    nb_size = 10, 
                    cohort = NULL, 
                    pct_drop = 1/10000,   
                    min_prob_increase = 0.05, max_iters = 40, logresults = FALSE) {

  #### preliminaries -----------------------------------
  # infer bg if not provided: assume background is proportional to the scaling factor s
  s <- rowSums(counts)
  if (is.null(bg)) {
    bgmod <- stats::lm(neg ~ s - 1)
    bg <- bgmod$fitted
  }
  if (length(bg) == 1) {
    bg = rep(bg, nrow(counts))
  }
  if (!is.null(fixed_profiles)) {
    if (!identical(rownames(fixed_profiles), rownames(counts))) {
      stop("gene ids in fixed profiles and counts aren't aligned")
    }
  }

  clusterlog = NULL

  if (is.null(cohort)) {
    cohort <- rep("all", length(neg))
  }
  
  #### get initial profiles: ----------------------------------
  
  if (is.null(init_profiles) & is.null(init_clust)) {
    stop("Must specify either init_clust or init_profiles")
  }
    

  # if init_profiles are provided, take them:
  if (!is.null(init_profiles)) {
    clust_old <- rep("unassigned", nrow(counts))
    names(clust_old) <- rownames(counts)
    profiles <- init_profiles
    clustnames <- colnames(init_profiles)
  } 
  # if no init_profiles are provided, derive them:
  if (is.null(init_profiles)) {
    clust_old = init_clust
    names(clust_old) <- rownames(counts)
    # derive first profiles from init_clust
    profiles <- Estep(counts = counts[!is.na(clust_old), ],
                      clust = init_clust[!is.na(clust_old)],
                      neg = neg[!is.na(clust_old)])
    # keep fixed_profiles unchanged:
    profiles[, colnames(fixed_profiles), drop = FALSE] <- fixed_profiles
    clustnames <- colnames(profiles)
  }
  
  #### run EM algorithm iterations: ----------------------------------
  pct_changed = c()
  if (logresults) {
    clusterlog <- init_clust
  }
  
  for (iter in seq_len(max_iters)) {
    message(paste0("iter ", iter))
    # M-step: get cell * cluster probs:
    
    probs <- Mstep(counts = counts,
                   means = profiles,
                   cohort = cohort, 
                   bg = bg,
                   size = nb_size)
    if (logresults) {
      clusterlog <- cbind(clusterlog, colnames(probs)[apply(probs, 1, which.max)])
    }
    
    oldprofiles <- profiles

    # E-step: update profiles:
    if (n_clusts != 0){
      tempclust <- colnames(probs)[apply(probs, 1, which.max)]
      profiles <- Estep(counts = counts,
                        clust = tempclust,
                        neg = neg)        
    }
    
    # for any profiles that have been lost, replace them with their previous version:
    lostprofiles = names(which(colSums(!is.na(profiles)) == 0))
    profiles[, lostprofiles] = oldprofiles[, lostprofiles]
    # keep fixed_profiles unchanged:
    profiles[, colnames(fixed_profiles), drop = FALSE] <- fixed_profiles
    
    # get cluster assignment
    clust = colnames(probs)[apply(probs, 1, which.max)]
    if (iter == 1){
      pct_changed <- mean(clust != clust_old)
    } else {
      index_changes <- which(clust != clust_old)
      probs_max <- apply(probs, 1, max)
      index_valid_changes <- which((probs_max - probs_old_max)[index_changes] > min_prob_increase)

      # record number of switches
      pct_changed = c(pct_changed, round(length(index_valid_changes)/nrow(probs), 5))

      if ( length(index_valid_changes)/nrow(probs) <= pct_drop ){
        message(sprintf("Converged: <= %s%% of cell type assignments changed in the last iteration.", pct_drop*100))
        message(        "==========================================================================")
        break
      }
    }
    clust_old = colnames(probs)[apply(probs, 1, which.max)]
    probs_old_max = apply(probs, 1, max)
  }
  names(pct_changed) <- paste0("Iter_", seq_len(iter))

  out = list(clust = clust,
             probs = probs,
             profiles = sweep(profiles, 2, colSums(profiles), "/") * 1000,
             pct_changed = pct_changed,
             #logliks = logliks,
             clusterlog = clusterlog)
  return(out)
}

#' For a numeric object, return a logical object of whether each element is the max or not.
#' @param x a vector of values
#' @return a vecetor of logical values
#'
ismax <- function(x) {
  return(x == max(x, na.rm = T))
}



#' Update logliks based on frequencies
#' @param logliks Matrix of cells' (rows) loglikelihoods under clusters (columns)
#' @param cohort Vector of cells' cohort memberships
#' @param minfreq Minimum frequency to give any cell type in any cohort
#' @param nbaselinecells Number of cells from baseline distribution to add to the 
#'  cohort-specific frequencies, thereby shrinking each cohort's data towards the population
#' @return An adjusted logliks matrix
update_logliks_with_cohort_freqs <- function(logliks, cohort, minfreq = 1e-6, nbaselinecells = 100) {
  # get overall cluster frequencies:
  clust <- colnames(logliks)[apply(logliks, 1, which.max)]
  baselinefreqs <- prop.table(table(clust))
  baselinefreqs[setdiff(unique(colnames(logliks)), names(baselinefreqs))] <- 0
  baselinefreqs <- baselinefreqs[colnames(logliks)]
  
  for (cohortname in unique(cohort)) {
    use <- cohort == cohortname
    # get cluster frequencies in cohort:
    cohortabundances <- (table(clust[use]))
    cohortabundances[setdiff(unique(colnames(logliks)), names(cohortabundances))] <- 0
    cohortabundances <- cohortabundances[colnames(logliks)]
    
    # add atop 100 cells' worth of baselinefreqs:
    cohortfreqs <- prop.table(baselinefreqs * nbaselinecells + cohortabundances)
    cohortfreqs <- cohortfreqs + minfreq
    # adjust the logliks:
    logliks[use, ] <- sweep(logliks[use, , drop = FALSE], 2, log(cohortfreqs), "+")
  }
  return(logliks)
}
