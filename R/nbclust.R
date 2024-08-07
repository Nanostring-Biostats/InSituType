#' Calculate the likelihood of the expression mat
#'   using the reference profiles of x
#'   
#' @param x a vector of a reference mean profile for the cell type
#' @param xsd a vector of a reference standard deviation profile for the cell type
#' @param mat a matrix of expression levels in all cells: for Protein data, we use raw data for calculating the scaling factor 
#' @param bg background level (default: 0.01)
#' @param size the parameters for dnbinom function (default: 10)
#' @param digits the number of digits for rounding
#' @param assay_type Assay type of RNA, protein (default = "rna")
#'
#' @importFrom Matrix rowSums
#' @importFrom stats dnbinom
#' @importFrom methods as is
#' @return likelihood for profile
#' @return cells x profiles matrix of log likelihoods
#' @examples
#' data("mini_nsclc")
#' data("ioprofiles")
#' bg <- Matrix::rowMeans(mini_nsclc$neg)
#' genes <- intersect(dimnames(mini_nsclc$counts)[[2]], dimnames(ioprofiles)[[1]])
#' mat <- mini_nsclc$counts[, genes]
#' x <- ioprofiles[genes, ]
#' lldist(x = x, mat = mini_nsclc$counts, bg = bg, assay_type="RNA")
#' 

lldist <- function(x, xsd=NULL, mat, 
                   bg = 0.01, size = 10, digits = 2, 
                   assay_type = c("rna", "protein")) {
  
  assay_type <- match.arg(tolower(assay_type), c("rna", "protein"))

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
  # Check dimensions on bg and stop with informative error if not
  #  conformant
  if (is.vector(bg)) {
    if (!identical(length(bg), nrow(mat))) {
      errorMessage <- sprintf("Dimensions of count matrix and background are not conformant.\nCount matrix rows: %d, length of bg: %d",
                              nrow(mat), length(bg))
      stop(errorMessage)
    }
    
    bgsub <- mat
    bgsub@x <- bgsub@x - bg[bgsub@i + 1]
    bgsub@x <- pmax(bgsub@x, 0)
    bgsub <- Matrix::rowSums(bgsub)
    
    if(identical(tolower(assay_type), "rna")){
      xsd <- NULL
      res <- lls_rna(mat=mat, bgsub=bgsub, x=x, bg=bg, size=size)
    }
    
    if(identical(tolower(assay_type), "protein")){
      res <- lls_protein(mat=as.matrix(mat), bgsub=bgsub, x=as.matrix(x), xsd=as.matrix(xsd))
    }
    
  }else{
    # non-optimized code used if bg is cell x gene matrix
    bgsub <- pmax(mat - bg, 0)
    
    sum_of_x <- sum(x)
    s <- Matrix::rowSums(bgsub) / sum_of_x
    
    # override it if s is negative:
    s[s <= 0] <- Matrix::rowSums(mat[s <= 0, , drop = FALSE]) / sum_of_x
    
    if(identical(tolower(assay_type), "rna")){
      yhat <- s %*% t(x) + bg
      res <- stats::dnbinom(x = as.matrix(mat), size = size, mu = yhat, log = TRUE)
    }
    
    if(identical(tolower(assay_type), "protein")){
      yhat <- s %*% t(x)
      ysd <- s %*% t(xsd)
      ## Estimate SD for each Cell type and each Protein. 
      ## Use the SD for all cells each protein with the same cell type although yhat is different across cells. means are more shift, but sd is overall distribution shape
      
      res <- stats::dnorm(x = as.matrix(mat), sd = ysd, mean = yhat, log = TRUE)
      if(is.matrix(res)){
        res <- Matrix::rowSums(res)
      }
    }
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
#' @param sds Matrix of standard deviation cluster profiles,
#'  with genes in rows and clusters in columns.
#' @param cohort a vector of cells' "cohort" assignment, used to update logliks 
#'  based on cluster frequencies within a cohort.
#' @param bg Expected background
#' @param size NB size parameter
#' @param digits Round the output to this many digits (saves memory)
#' @param return_loglik If TRUE, logliks will be returned. If FALSE, probabilities will be returned. 
#' @param assay_type Assay type of RNA, protein (default = "rna")
#' 
#' @return Matrix of probabilities of each cell belonging to each cluster
#' @export
#' @examples 
#' data("mini_nsclc")
#' data("ioprofiles")
#' sharedgenes <- intersect(rownames(ioprofiles), colnames(mini_nsclc$counts))
#' Mstep(mini_nsclc$counts, ioprofiles[sharedgenes, ], bg = Matrix::rowMeans(mini_nsclc$neg), cohort = NULL, assay_type="RNA")
  
Mstep <- function(counts, means, sds=NULL, 
                  cohort, bg = 0.01, size = 10, 
                  digits = 2, return_loglik = FALSE, 
                  assay_type = c("rna", "protein")) {
  assay_type <- match.arg(tolower(assay_type), c("rna", "protein"))
  
  # get logliks of cells * clusters
  logliks <- lldist(x = means,
                    mat = counts,
                    xsd = sds,
                    bg = bg,
                    size = size,
                    digits = digits,
                    assay_type = assay_type)
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
#' @param assay_type Assay type of RNA, protein (default = "rna")
#'
#' @importFrom Matrix rowSums
#'
#' @return A list with two elements: 1.  A matrix of cluster profiles, genes * clusters. 
#'         2. In protein mode, a matrix holding SDs, also genes * clusters. NULL in RNA mode.
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
#'  max_iters = 5,
#'  assay_type="RNA",
#' ) # choosing inadvisably low numbers to speed the vignette; using the defaults in recommended.
#' Estep(counts = mini_nsclc$counts, clust = unsup$clust, neg = Matrix::rowMeans(mini_nsclc$neg), assay_type="RNA")

Estep <- function(counts, clust, neg, 
                  assay_type = c("rna", "protein")) {
  
  assay_type <- match.arg(tolower(assay_type), c("rna", "protein"))

  # get cluster means:
  means <- sapply(unique(clust), function(cl) {
    
    if(identical(tolower(assay_type), "rna")){
      means = pmax(Matrix::colMeans(counts[clust == cl, , drop = FALSE]) - mean(neg[clust == cl]), 0)
    }
    
    if(identical(tolower(assay_type), "protein")){
      means = Matrix::colMeans(counts[clust == cl, , drop = FALSE])  #- mean(neg[clust == cl])
    }
    return(means)
  })

  if(identical(tolower(assay_type), "rna")){
    sds = NULL
  }
  if(identical(tolower(assay_type), "protein")){
    sds <- sapply(unique(clust), function(cl) {
      sds = apply(counts[clust == cl, , drop = FALSE], 2, sd)  #- sd(neg[clust == cl])
      return(sds)
    })
  }
  if (is.matrix(sds)) {
    rownames(sds) <- rownames(means)
  }
  
  return(list(profiles=means, sds=sds))
}


#' Cluster via EM algorithm based on cell logliks
#'
#' Cluster single cell gene expression data using an EM algorithm.
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negative probe counts per cell. 
#' @param assay_type Assay type of RNA, protein (default = "rna")
#' @param bg Expected background
#' @param fixed_profiles Matrix of mean expression profiles to hold unchanged throughout iterations. genes * cell types
#' @param fixed_sds Matrix of standard deviation profiles of pre-defined
#'   clusters to hold unchanged throughout iterations. 
#'   Columns must all be included in the init_clust variable. This parameter is 
#'   only for assay_type of protein.
#' @param init_profiles Matrix of cluster mean profiles under which to begin iterations.
#' If NULL, initial assignments will be automatically inferred, using init_clust 
#' if available, and using random clusters if not. 
#' @param init_sds Matrix of cluster SDs profiles under which to begin iterations.
#' If NULL, initial assignments will be automatically inferred, using init_clust 
#' if available, and using random clusters if not. Only for assay_type of protein
#' @param init_clust Vector of initial cluster assignments.
#' If NULL, initial assignments will be automatically inferred.
#' @param nb_size The size parameter to assume for the NB distribution. Only for assay_type of RNA.
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
#'        neg =  Matrix::rowMeans(mini_nsclc$neg), 
#'        assay_type = "RNA",
#'        bg = NULL,
#'        fixed_profiles = ioprofiles[sharedgenes, 1:3],
#'        init_profiles = NULL, 
#'        init_clust = rep(c("a", "b"), 
#'        nrow(mini_nsclc$counts) / 2),
#'        nb_size = 10,
#'        cohort = rep("a", nrow(mini_nsclc$counts)),
#'        pct_drop = 1/10000,
#'        min_prob_increase = 0.05, 
#'        max_iters = 3, 
#'        logresults = FALSE)
nbclust <- function(counts, 
                    neg, 
                    assay_type = c("rna", "protein"), 
                    bg = NULL, 
                    fixed_profiles = NULL,
                    fixed_sds = NULL,
                    init_profiles = NULL, 
                    init_sds = NULL, 
                    init_clust = NULL, 
                    nb_size = 10,
                    cohort = NULL, 
                    pct_drop = 1/10000,   
                    min_prob_increase = 0.05, 
                    max_iters = 40, 
                    logresults = FALSE) {
  assay_type <- match.arg(tolower(assay_type), c("rna", "protein"))

  #### preliminaries -----------------------------------
  # infer bg if not provided: assume background is proportional to the scaling factor s
  
  # get vector of expected background:
  bg <- estimateBackground(counts = counts, neg = neg, bg = bg)
  

  if (!is.null(fixed_profiles)) {
    if (!identical(rownames(fixed_profiles), colnames(counts))) {
      stop("gene ids in fixed profiles and counts aren't aligned")
    }
  }

  clusterlog <- NULL

  if (is.null(cohort)) {
    cohort <- rep("all", length(bg))
  }
  
  #### get initial profiles: ----------------------------------
  
  if ((is.null(init_profiles) && is.null(init_clust)) || (is.null(init_profiles) && is.null(init_sds) && is.null(init_clust)) ) {
    stop("Must specify either init_clust or init_profiles")
  }

  # if init_profiles are provided, take them:
  if (!is.null(init_profiles)) {
    clust_old <- rep("unassigned", nrow(counts))
    names(clust_old) <- rownames(counts)
    profiles <- init_profiles
    sds <- init_sds
  } 
  # if no init_profiles for RNA & protein, derive them:
  if (is.null(init_profiles)) {
    clust_old <- init_clust
    names(clust_old) <- rownames(counts)
    # derive first profiles from init_clust
    profiles_info <- Estep(counts = counts[!is.na(clust_old), ],
                      clust = init_clust[!is.na(clust_old)],
                      neg = bg[!is.na(clust_old)],
                      assay_type=assay_type)
    ## if no init_profiles is provided, use the derived profiles
    profiles <- profiles_info$profiles
    sds <- profiles_info$sds
    
    ## if no init_sds is provided, use the derived SD profiles
    ## otherwise, use the init_sds
    if(identical(tolower(assay_type), "rna")){
      sds <- NULL
    }
  }
  
  # keep fixed_profiles unchanged:
  if (length(profiles) == 0) {
    profiles <- NULL
  }
  profiles <- cbind(profiles[, setdiff(colnames(profiles), colnames(fixed_profiles)), drop = FALSE], fixed_profiles)
  
  if(identical(tolower(assay_type), "protein")){
    sds <- cbind(sds[, setdiff(colnames(sds), colnames(fixed_sds)), drop = FALSE], fixed_sds)
  }
  clustnames <- colnames(profiles)

  #### run EM algorithm iterations: ----------------------------------
  pct_changed <- c()
  if (logresults) {
    clusterlog <- init_clust
  }
  
  for (iter in seq_len(max_iters)) {
    if (iter %% 5 == 0) {
      message(paste0("iter ", iter))
    }
    # M-step: get cell * cluster probs:
    probs <- Mstep(counts = counts,
                   means = profiles,
                   sds=sds,
                   cohort = cohort, 
                   bg = bg,
                   size = nb_size,
                   assay_type=assay_type)
    if (logresults) {
      clusterlog <- cbind(clusterlog, colnames(probs)[apply(probs, 1, which.max)])
    }
    
    oldprofiles <- profiles
    oldsds <- sds
    
    # E-step: update profiles:
    tempclust <- colnames(probs)[apply(probs, 1, which.max)]
    profiles_info <- Estep(counts = counts,
                      clust = tempclust,
                      neg = bg, 
                      assay_type=assay_type)
    
    profiles <- profiles_info$profiles
    sds <- profiles_info$sds
    
    # for any profiles that have been lost, replace them with their previous version:
    lostprofiles <- setdiff(clustnames, colnames(profiles))
    profiles <- cbind(profiles, oldprofiles[, lostprofiles, drop = FALSE])[, clustnames]
    
    if(identical(tolower(assay_type), "rna")){
      sds <- NULL
    }
    if(identical(tolower(assay_type), "protein")){
      sds <- cbind(sds, oldsds[, lostprofiles, drop = FALSE])[, clustnames]
    }
    
    # keep fixed_profiles unchanged:
    profiles[, colnames(fixed_profiles)] <- as.vector(fixed_profiles)
    
    if(identical(tolower(assay_type), "rna")){
      sds <- NULL
    }
    if(identical(tolower(assay_type), "protein")){
      sds[, colnames(fixed_sds)] <- as.vector(fixed_sds)
    }
    
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

  if(identical(tolower(assay_type), "protein")){
    sds = sweep(sds, 2, colSums(sds), "/") * 1000
  }
  
  out <- list(clust = clust,
             probs = probs,
             profiles = sweep(profiles, 2, colSums(profiles), "/") * 1000,
             sds = sds,
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
