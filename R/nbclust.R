#' Calculate the likelihood of the expression mat
#'   using the reference profiles of x
#'
#' @param x matrix of reference cell types, genes x profiles
#' @param mat a matrix of expression levels in all cells, cells x genes
#' @param bg background level (default: 0.01), usually vector length equal to
#'   number of cells
#' @param size the parameters for dnbinom function (default: 10)
#' @param digits the number of digits for rounding
#'
#' @importFrom Matrix rowSums
#' @importFrom stats dnbinom
#' @importFrom methods as is
#' @return cells x profiles matrix of log likelihoods
#' @examples
#' data("mini_nsclc")
#' data("ioprofiles")
#' bg <- Matrix::rowMeans(mini_nsclc$neg)
#' genes <- intersect(dimnames(mini_nsclc$counts)[[2]], dimnames(ioprofiles)[[1]])
#' mat <- mini_nsclc$counts[, genes]
#' x <- ioprofiles[genes, ]
#' lldist(x = x, mat = mini_nsclc$counts, bg = Matrix::rowMeans(mini_nsclc$neg))
lldist <- function(x, mat, bg = 0.01, size = 10, digits = 2) {
  # convert to matrix form if only a vector was input:
  if (is.vector(mat)) {
    mat <- as(matrix(mat, nrow = 1), "dgCMatrix")
  } else if (is.matrix(mat)) {
    mat <- as(mat, "dgCMatrix")
  } else if (!is(mat, "dgCMatrix")) {
    errorMessage <-
      sprintf(
        "The `type` of parameter `mat` needs to be of one of dgCMatrix, vector, matrix, array, but is found to be of type %s",
        class(mat)
      )
    stop(errorMessage)
  }
  if (is.vector(bg)) {
    # Check dimensions on bg and stop with informative error if not
    #  conformant
    if (!identical(length(bg), nrow(mat))) {
      errorMessage <- sprintf("Dimensions of count matrix and background are not conformant.\nCount matrix rows: %d, length of bg: %d",
                              nrow(mat), length(bg))
      stop(errorMessage)
    }
    bgsub <- mat
    bgsub@x <- bgsub@x - bg[bgsub@i + 1]
    bgsub@x <- pmax(bgsub@x, 0)
    bgsub <- Matrix::rowSums(bgsub)
    res <- fast_lldist(mat, bgsub, x, bg, size)
  } else {
    # non-optimized code used if bg is cell x gene matrix
    bgsub <- pmax(mat - bg, 0)
    res <- apply(x, 2, function(profile) {
      sum_of_x <- sum(profile)
      s <- Matrix::rowSums(bgsub) / sum_of_x
      # override it if s is negative:
      s[s <= 0] <- Matrix::rowSums(mat[s <= 0, , drop = FALSE]) / sum_of_x
      yhat <- s %*% t(profile) + bg
      rowSums(stats::dnbinom(x = as.matrix(mat), size = size, mu = yhat, log = TRUE))
    })
  }

  rownames(res) <- rownames(mat)
  colnames(res) <- colnames(x)
  return(round(res, digits))
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
#' @examples 
#' data("mini_nsclc")
#' data("ioprofiles")
#' sharedgenes <- intersect(rownames(ioprofiles), colnames(mini_nsclc$counts))
#' Mstep(mini_nsclc$counts, ioprofiles[sharedgenes, ], bg = Matrix::rowMeans(mini_nsclc$neg), cohort = NULL)
Mstep <- function(counts, means, cohort, bg = 0.01, size = 10, digits = 2, return_loglik = FALSE) {
  # get logliks of cells * clusters
  logliks <- lldist(x = means,
                    mat = counts,
                    bg = bg,
                    size = size,
                    digits = digits)
  
  # adjust by cohort frequency:
  logliks <- update_logliks_with_cohort_freqs(logliks = logliks, 
                                              cohort = cohort, 
                                              minfreq = 1e-4, 
                                              nbaselinecells = 100) 
  if (return_loglik) {
    return(round(logliks, digits))
  } else {
    # first rescale (ie recenter on log scale) to avoid rounding errors:
    logliks <- sweep(logliks, 1, apply(logliks, 1, max), "-")
    # get on likelihood scale:
    liks <- exp(logliks)
    # convert to probs
    probs <- sweep(liks, 1, rowSums(liks), "/")
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
#' @examples 
#' data("ioprofiles")
#' unsup <- insitutype(
#'  x = mini_nsclc$counts,
#'  neg = Matrix::rowMeans(mini_nsclc$neg),
#'  n_clusts = 8,
#'  n_phase1 = 200,
#'  n_phase2 = 500,
#'  n_phase3 = 2000,
#'  n_starts = 1,
#'  max_iters = 5
#' ) # choosing inadvisably low numbers to speed the vignette; using the defaults in recommended.
#' Estep(counts = mini_nsclc$counts, clust = unsup$clust, neg = Matrix::rowMeans(mini_nsclc$neg))
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
#' @param fixed_profiles Matrix of cluster profiles to hold unchanged throughout iterations.
#' @param init_profiles Matrix of cluster profiles under which to begin iterations.
#' If NULL, initial assignments will be automatically inferred, using init_clust 
#' if available, and using random clusters if not. 
#' @param init_clust Vector of initial cluster assignments.
#' If NULL, initial assignments will be automatically inferred.
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param cohort Vector of cells' "cohort" assignments, uses to assess frequencies in each cluster. 
#' @param pct_drop the decrease in percentage of cell types with a valid switchover to 
#'  another cell type compared to the last iteration. Default value: 1/10000. A valid 
#'  switchover is only applicable when a cell has changed the assigned cell type with its
#'  highest cell type probability increased by min_prob_increase. 
#' @param min_prob_increase the threshold of probability used to determine a valid cell 
#'  type switchover
#' @param max_iters Maximum number of iterations
#' @param logresults Populate clusterlog in returned list
#'
#'  @importFrom stats lm
#'
#' @return A list, with the following elements:
#' \enumerate{
#' \item probs: a matrix of probabilities of all cells (rows) belonging to all clusters (columns)
#' \item profiles: a matrix of cluster-specific expression profiles
#' }
#' @examples 
#' data("ioprofiles")
#' data("mini_nsclc")
#' sharedgenes <- intersect(colnames(mini_nsclc$counts), rownames(ioprofiles))
#' nbclust(counts = mini_nsclc$counts[, sharedgenes],
#'        neg =  Matrix::rowMeans(mini_nsclc$neg), bg = NULL,
#'        fixed_profiles = ioprofiles[sharedgenes, 1:3],
#'        init_profiles = NULL, init_clust = rep(c("a", "b"), nrow(mini_nsclc$counts) / 2),
#'        nb_size = 10,
#'        cohort = rep("a", nrow(mini_nsclc$counts)),
#'        pct_drop = 1/10000,
#'        min_prob_increase = 0.05, max_iters = 3, logresults = FALSE)
nbclust <- function(counts, neg, bg = NULL, 
                    fixed_profiles = NULL,
                    init_profiles = NULL, init_clust = NULL, 
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
    bg <- rep(bg, nrow(counts))
  }
  if (!is.null(fixed_profiles)) {
    if (!identical(rownames(fixed_profiles), colnames(counts))) {
      stop("gene ids in fixed profiles and counts aren't aligned")
    }
  }

  clusterlog <- NULL

  if (is.null(cohort)) {
    cohort <- rep("all", length(neg))
  }
  
  #### get initial profiles: ----------------------------------
  
  if (is.null(init_profiles) && is.null(init_clust)) {
    stop("Must specify either init_clust or init_profiles")
  }

  # if init_profiles are provided, take them:
  if (!is.null(init_profiles)) {
    clust_old <- rep("unassigned", nrow(counts))
    names(clust_old) <- rownames(counts)
    profiles <- init_profiles
  } 
  # if no init_profiles are provided, derive them:
  if (is.null(init_profiles)) {
    clust_old <- init_clust
    names(clust_old) <- rownames(counts)
    # derive first profiles from init_clust
    profiles <- Estep(counts = counts[!is.na(clust_old), ],
                      clust = init_clust[!is.na(clust_old)],
                      neg = neg[!is.na(clust_old)])
  }
  # keep fixed_profiles unchanged:
  if (length(profiles) == 0) {
    profiles <- NULL
  }
  profiles <- cbind(profiles[, setdiff(colnames(profiles), colnames(fixed_profiles)), drop = FALSE], fixed_profiles)
  clustnames <- colnames(profiles)
  
  #### run EM algorithm iterations: ----------------------------------
  pct_changed <- c()
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
    tempclust <- colnames(probs)[apply(probs, 1, which.max)]
    profiles <- Estep(counts = counts,
                      clust = tempclust,
                      neg = neg)        
    
    # for any profiles that have been lost, replace them with their previous version:
    lostprofiles <- setdiff(clustnames, colnames(profiles))
    profiles <- cbind(profiles, oldprofiles[, lostprofiles, drop = FALSE])[, clustnames]
    # keep fixed_profiles unchanged:
    profiles[, colnames(fixed_profiles)] <- as.vector(fixed_profiles)
    
    # get cluster assignment
    clust <- colnames(probs)[apply(probs, 1, which.max)]
    if (iter == 1) {
      pct_changed <- mean(clust != clust_old)
    } else {
      index_changes <- which(clust != clust_old)
      probs_max <- apply(probs, 1, max)
      index_valid_changes <- which((probs_max - probs_old_max)[index_changes] > min_prob_increase)

      # record number of switches
      pct_changed <- c(pct_changed, round(length(index_valid_changes) / nrow(probs), 5))

      if (length(index_valid_changes) / nrow(probs) <= pct_drop) {
        message(sprintf("Converged: <= %s%% of cell type assignments changed in the last iteration.", pct_drop * 100))
        message("==========================================================================")
        break
      }
    }
    clust_old <- colnames(probs)[apply(probs, 1, which.max)]
    probs_old_max <- apply(probs, 1, max)
  }
  names(pct_changed) <- paste0("Iter_", seq_len(iter))

  out <- list(clust = clust,
             probs = probs,
             profiles = sweep(profiles, 2, colSums(profiles), "/") * 1000,
             pct_changed = pct_changed,
             clusterlog = clusterlog)
  return(out)
}

#' For a numeric object, return a logical object of whether each element is the max or not.
#' @param x a vector of values
#' @return a vecetor of logical values
#' @examples
#' ismax(c(3, 5, 5, 2))
ismax <- function(x) {
  return(x == max(x, na.rm = TRUE))
}



#' Update logliks based on frequencies
#' @param logliks Matrix of cells' (rows) loglikelihoods under clusters (columns)
#' @param cohort Vector of cells' cohort memberships
#' @param minfreq Minimum frequency to give any cell type in any cohort
#' @param nbaselinecells Number of cells from baseline distribution to add to the 
#'  cohort-specific frequencies, thereby shrinking each cohort's data towards the population
#' @return An adjusted logliks matrix
update_logliks_with_cohort_freqs <- function(logliks, cohort, minfreq = 1e-6, nbaselinecells = 50) {
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
