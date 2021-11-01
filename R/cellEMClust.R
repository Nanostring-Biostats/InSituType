
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
#' @param bg Expected background
#' @param size NB size parameter
#' @param digits Round the output to this many digits (saves memory)
#' @return Matrix of probabilities of each cell belonging to each cluster
Mstep <- function(counts, means, bg = 0.01, size = 10, digits = 2) {
  # get logliks of cells * clusters
  logliks <- apply(means, 2, function(x) {
    lldist(x = x, mat = counts, bg = bg, size = size)
  })
  # first rescale (ie recenter on log scale) to avoid rounding errors:
  logliks <- sweep( logliks , 1 , apply( logliks , 1 , max ) , "-" )
  # get on likelihood scale:
  liks <- exp( logliks )
  # convert to probs
  probs <- sweep( liks , 1 , rowSums( liks ) , "/" )

  return(round(probs, digits))
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
#'
Estep <- function(counts, clust, neg) {

  # get cluster means:
  if (is.vector(clust)) {
    means <- sapply(unique(clust), function(cl) {
      pmax(Matrix::colSums(counts[clust == cl, , drop = FALSE]) - sum(neg[clust == cl]), 0)
    })
  }
  if (is.matrix(clust)) {
    means <- apply(clust, 2, function(x) {
      wts <- x / sum(x)
      pmax(Matrix::colSums(sweep(counts, 1, wts, "*")) - sum(neg * wts), 0)
    })
  }

  return(means)
}


#' E step: update reference cell profiles
#'
#' Given cell assignments (or posterior probabilities), update the pre-defined
#'  reference cell profiles.
#' @param counts Counts matrix, cells * genes.
#' @param clust Vector of cluster assignments, or a matrix of probabilities
#'   of cells (rows) belonging to clusters (columns).
#' @param neg Vector of mean background counts
#' @param fixed_profiles Matrix of expression profiles of pre-defined clusters,
#'  e.g. from previous scRNA-seq. These profiles will not be updated by the EM algorithm.
#'  Colnames must all be included in the init_clust variable.
#'
#' @importFrom Matrix rowSums
#'
#' @return A matrix of cluster profiles, genes * clusters
#'
Estep_reference <- function(counts, clust, neg, fixed_profiles) {


  # get cluster means:
  # if clust has been passed as a vector of cluster names:
  if (is.vector(clust)) {
    means = sapply(unique(clust), function(cl) {
      pmax(Matrix::colSums(counts[clust == cl, , drop = FALSE]) - sum(neg[clust == cl]), 0)
    })
  }
  # if clust has been passed as a matrix of cluster probabilities:
  if (is.matrix(clust)) {
    means <- apply(clust, 2, function(x) {
      wts <- x / sum(x)
      pmax(Matrix::colSums(sweep(counts, 1, wts, "*")) - sum(neg * wts), 0)
    })
  }

  # estimate gene scaling effects
  gene_scaling_factors <- estimate_platform_effects(
    mean_profiles = means,
    fixed_profiles = fixed_profiles
  )

  # rescale fixed profiles by platform effects:
  updated_profiles <- sweep(fixed_profiles, 1, gene_scaling_factors, "*")

  return(updated_profiles)
}


#' Estimate NB size parameter
#'
#' Given cell assignments (or posterior probabilities), update theta.
#' Assumption is one theta for all data.
#' @param counts Counts matrix, cells * genes.
#' @param clust Vector of cluster assignments, or a matrix of probabilities
#'   of cells (rows) belonging to clusters (columns).
#' @param bg Expected background
#' @return A scalar giving the mle for the size parameter
Estep_size <- function(counts, clust, bg) {

  # define the matrix of expected counts (mu in the NB model):
  if (is.vector(clust)) {
    expected = NOTIMPLEMENTEDYET

  }
  if (is.matrix(clust)) {
    expected = clust %*% means + bg
    # (nope, have to treat bg mroe carefully. might be vec or matrix)
    # also needs to incorporate s somehow
  }

  # fit theta/size:
  mod = glm.nb(as.vector(counts) ~ as.vector(expected))$theta
  return(theta)
}

#' Cluster via EM algorithm based on cell logliks
#'
#' Cluster single cell gene expression data using an EM algorithm.
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param bg Expected background
#' @param fixed_profiles Matrix of expression profiles of pre-defined clusters,
#'  e.g. from previous scRNA-seq. These profiles will not be updated by the EM algorithm.
#'  Colnames must all be included in the init_clust variable.
#' @param init_free_profiles Matrix of cluster profiles under which to begin iterations.
#' Should NOT contain the fixed_profiles.
#' If NULL, initial assignments will be automatically inferred, using init_clust 
#' if available, and using random clusters if not. 
#' @param init_clust Vector of initial cluster assignments.
#' If NULL, initial assignments will be automatically inferred.
#' @param n_clusts Number of clusters, in addition to any pre-specified cell types.
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param n_iters Numer of iterations
#' @param method Whether to run a SEM algorithm (points given a probability
#'   of being in each cluster) or a classic EM algorithm (all points wholly
#'   assigned to one cluster).
#' @param shrinkage Fraction by which to shrink the average profiles towards
#'  the fixed profiles. 1 = keep the fixed profile; 0 = don't shrink the mean profile.
#' @param updated_reference Rescaled, possibly shrunken version of fixed_profiles.
#'  This argument is intended to be used by cellEMclust, not by the user.
#' @param pct_drop the decrease in percentage of cell types with a valid switchover to 
#'  another cell type compared to the last iteration. Default value: 1/10000. A valid 
#'  switchover is only applicable when a cell has changed the assigned cell type with its
#'  highest cell type probability increased by min_prob_increase. 
#' @param min_prob_increase the threshold of probability used to determine a valid cell 
#'  type switchover
#'
#'  @importFrom stats lm
#'
#' @return A list, with the following elements:
#' \enumerate{
#' \item probs: a matrix of probabilities of all cells (rows) belonging to all clusters (columns)
#' \item profiles: a matrix of cluster-specific expression profiles
#' }
nbclust <- function(counts, neg, bg = NULL, 
                    init_free_profiles = NULL, init_clust = NULL, n_clusts = NULL,
                    fixed_profiles = NULL, nb_size = 10, n_iters = 20,
                    method = "CEM", shrinkage = 0.8,
                    updated_reference = NULL, pct_drop = 1/10000,   
                    min_prob_increase = 0.05) {

  #### preliminaries -----------------------------------
  # infer bg if not provided: assume background is proportional to the scaling factor s
  if (is.null(bg)) {
    s <- rowSums(counts)
    bgmod <- stats::lm(neg ~ s - 1)
    bg <- bgmod$fitted
  }

  # specify which profiles to leave fixed:
  if (!is.null(fixed_profiles)) {
    fixed_profile_names = colnames(fixed_profiles)
  } else {
    fixed_profile_names <- NULL
  }

  # start with updated_reference = fixed_profiles if not already specified
  if (is.null(updated_reference)) {    
    updated_reference <- fixed_profiles
  }


  #### get initial profiles: ----------------------------------

  # if init_free_profiles are provided, take them:
  if (!is.null(init_free_profiles)) {
    free_profiles <- init_free_profiles
    clust_old <- rep("unassigned", nrow(counts))
    names(clust_old) <- rownames(counts)
  } 
  # if no init_free_profiles are provided, derive them:
  if (is.null(init_free_profiles)) {
    # if no initial clustering is available, define it randomly:
    if (is.null(init_clust)) {
      if (is.null(n_clusts)) {
        stop("Must specify either init_clust or n_clusts.")
      }
      if (!is.null(fixed_profiles)) {
        n_fixed_profiles <- ncol(fixed_profiles)
      } else {
        n_fixed_profiles <- 0
      }
      clustnames <- makeClusterNames( fixed_profile_names , n_clusts + n_fixed_profiles )
      free_profile_names <- setdiff(clustnames, fixed_profile_names)
      
      # arbitrary but non-random initialization:
      init_clust = rep(clustnames, ceiling(nrow(counts) / length(clustnames)))[seq_len(nrow(counts))]
    }
  
    # subset on only the cells that aren't part of a pre-specified cluster:
    tempuse = !is.element(init_clust, fixed_profile_names)    #<------------------ randomizing in blocks might help discover cell types that are together in space. This is the most null randomization possible. 
    # use this subset to derive first free_profiles from init_clust
    free_profiles <- Estep(counts = counts[tempuse, ],
                           clust = init_clust[tempuse],
                           neg = neg[tempuse])
    
    clust_old = init_clust
  }
  # append free and fixed profiles:
  profiles <- cbind(updated_reference, free_profiles)
  #if (n_clusts > 0 ){
  #  profiles <- cbind(updated_reference, free_profiles)
  #} else {
  #  profiles <- updated_reference
  #}
  
  
  #### run EM algorithm iterations: ----------------------------------
  pct_changed = c()
  for (iter in seq_len(n_iters)) {
    message(paste0("iter ", iter))
    # M-step: get cell * cluster probs:
    probs <- Mstep(counts = counts,
                   means = profiles,
                   bg = bg,
                   size = nb_size)

    # E-step: update profiles:
    if (method == "CEM") {
      tempprobs = probs[, free_profile_names, drop = FALSE]
      # update the new cluster profiles:
      if (n_clusts != 0){
        free_profiles <- Estep(counts = counts,
                              clust = tempprobs,
                              neg = neg)
      }
      # update the reference profiles / "fixed_profiles"
      if( !is.null(updated_reference) ){
        updated_reference <-
          Estep_reference(counts = counts,
                          clust = probs[, fixed_profile_names],
                          neg = neg,
                          fixed_profiles = fixed_profiles)
      }
    }
    if (method == "EM") {
      # update the new cluster profiles:
      if (n_clusts != 0){
        tempprobs = 1 * t(apply(probs[, free_profile_names, drop = FALSE], 1, ismax))
        free_profiles <- Estep(counts = counts,
                              clust = tempprobs,
                              neg = neg)        
      }
      # update the reference profiles / "fixed_profiles"
      if( !is.null(updated_reference) ){
        updated_reference <-
          Estep_reference(counts = counts,
                          clust = probs[, fixed_profile_names],
                          neg = neg,
                          fixed_profiles = fixed_profiles)
      }

      # for any profiles that have been lost, replace them with their previous version:
      if ( n_clusts!=0 ){
        lostprofiles = names(which(colSums(!is.na(free_profiles)) == 0))
        free_profiles[, lostprofiles] = profiles[, lostprofiles]
      }
    }


    if ( n_clusts!=0 ){
      profiles <- cbind(updated_reference, free_profiles)
    } else {
      profiles <- updated_reference
    }
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
  # get loglik of each cell:
  logliks = unlist(sapply(seq_len(nrow(counts)), function(i) {   #<-------------------- should just assume vector background
    if (length(bg) == 1) {
      bgtemp = bg
    }
    if (length(bg) == nrow(counts)) {
      bgtemp = bg[i]
    }
    if (is.matrix(bg)) {
      bgtemp = bg[i, ]
    }
    lldist(mat = counts[i, ], x = profiles[, clust[i]], bg = bgtemp, size = nb_size)
  }))


  out = list(clust = clust,
             probs = probs,
             profiles = sweep(profiles, 2, colSums(profiles), "/") * 1000,
             pct_changed = pct_changed,
             logliks = logliks,
             updated_reference = updated_reference)
  return(out)
}

#' For a numeric object, return a logical object of whether each element is the max or not.
#' @param x a vector of values
#' @return a vecetor of logical values
#'

ismax <- function(x) {
  return(x == max(x, na.rm = T))
}


#' makeClusterNames
#' @param cNames existing cluster names
#' @param nClust number of clusters
#' @return a vector of complete cluster names
#'

makeClusterNames <- function( cNames , nClust )
{
  if ( is.null( cNames ) )
  {
    cNames <- letters[seq_len( nClust )]
    extraClustNames <- NULL
  }
  else if ( nClust > length( cNames ) )
  {
    nDiff <- nClust - length( cNames )
    if ( nDiff > 26 )
    {
      repLetters <- floor( nDiff / 26 )
      extraLetters <- nDiff %% 26
      extraClustNames <- paste0( rep( letters[seq_len( repLetters )] , each = 26 ) , letters )
      if ( extraLetters > 0 )
      {
        extraClustNames <- c( extraClustNames , paste0( letters[repLetters + 1] , letters[seq_len( extraLetters )] ) )
      }
    }
    else
    {
      extraClustNames <- letters[seq_len(nDiff)]
    }
  }
  else
  {
    if ( length( cNames ) < nClust )
    {
      stop( "nClust is smaller than nClust.  Do not know how to subset." )
    }
    extraClustNames <- NULL
  }
  clustnames <- c( cNames , extraClustNames )
  # make sure unique name when ignore case
  dup_flag <- duplicated( tolower( clustnames ) )
  if ( any( dup_flag ) )
  {
    # update on the duplicated name
    clustnames[dup_flag] <- make.unique( tolower( clustnames[dup_flag] ) , sep = '.' )
  }
  return( clustnames )
}


#' Clustering wrapper function
#'
#' A wrapper for nbclust, to manage subsampling and multiple random starts
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param bg Expected background
#' @param init_clust Vector of initial cluster assignments.
#' If NULL, initial assignments will be automatically inferred.
#' @param n_clusts Number of clusters, in addition to any pre-specified cell types.
#' @param fixed_profiles Matrix of expression profiles of pre-defined clusters,
#'  e.g. from previous scRNA-seq. These profiles will not be updated by the EM algorithm.
#'  Colnames must all be included in the init_clust variable.
#' @param align_genes Logical, for whether to align the counts matrix and the fixed_profiles by gene ID.
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param n_iters Numer of iterations
#' @param method Whether to run a SEM algorithm (points given a probability
#'   of being in each cluster) or a classic EM algorithm (all points wholly
#'   assigned to one cluster).
#' @param shrinkage Fraction by which to shrink the average profiles towards
#'  the fixed profiles. 1 = keep the fixed profile; 0 = don't shrink the mean profile.
#' @param subset_size To speed computations, each iteration will use a subset of only this many cells.
#'  (The final iteration runs on all cells.) Set to NULL to use all cells in every iter.
#'  (This option has not yet been enabled.)
#' @param n_starts the number of iterations
#' @param n_benchmark_cells the number of cells for benchmarking
#' @param n_final_iters the number of iterations on the final step
#' @param pct_drop the decrease in percentage of cell types with a valid switchover to 
#'  another cell type compared to the last iteration. Default value: 1/10000. A valid 
#'  switchover is only applicable when a cell has changed the assigned cell type with its
#'  highest cell type probability increased by min_prob_increase. 
#' @param min_prob_increase the threshold of probability used to determine a valid cell 
#'  type switchover
#'  
#' @importFrom stats lm
#' @importFrom Matrix rowMeans
#' @importFrom Matrix colSums
#'
#' @return A list, with the following elements:
#' \enumerate{
#' \item cluster: a vector given cells' cluster assignments
#' \item probs: a matrix of probabilies of all cells (rows) belonging to all clusters (columns)
#' \item profiles: a matrix of cluster-specific expression profiles
#' \item pct_changed: how many cells changed class at each step
#' \item logliks: a matrix of each cell's log-likelihood under each cluster
#' }
#' @export
#' @examples
#' # load data ("raw" and "cellannot"):
#' data(ioprofiles)
#' # predict per-cell bbackground:
#' bgmodel = lm(rowSums(raw[, grepl("NegPrb", colnames(raw))]) ~ rowSums(raw) - 1)$coef
#' bg.predicted = rowSums(raw) * bgmodel
#' # run unsupervised clustering with several random starts:
#' unsup <- cellEMClust(counts = raw,
#'                      bg = bg.predicted,
#'                      init_clust = NULL, n_clusts = 12,
#'                      fixed_profiles = NULL,
#'                      nb_size = 10,
#'                      n_iters = 10,  # this is not enough
#'                      method = "EM", shrinkage = 0.5,
#'                      subset_size = 1000,  # smaller than ideal
#'                      n_starts = 4,
#'                      n_benchmark_cells = 500,
#'                      n_final_iters = 5)   # this is not enough
#' # plot clusters:
#' scols = c(cellcols, brewer.pal(12, "Set3")[-2], brewer.pal(8, "Set2"))[1:length(unique(unsup$clust))]
#' names(scols) = unique(unsup$clust)
#' plot(cellannot$x, cellannot$y, pch = 16, col = scols[unsup$clust])
#'
#' # view immune oncology cell profiles (in ptolemy package data):
#' head(ioprofiles)
#' usegenes = intersect(rownames(ioprofiles), colnames(raw))
#' semi <- cellEMClust(counts = raw[, usegenes],
#'                     bg = bg.predicted,
#'                     init_clust = NULL, n_clusts = 6,
#'                     fixed_profiles = ioprofiles[usegenes, ],
#'                     nb_size = 10,
#'                     n_iters = 10,  # this is not enough
#'                     method = "EM",
#'                     shrinkage = 0.5,
#'                     subset_size = 1000,
#'                     n_starts = 4,
#'                     n_benchmark_cells = 500,
#'                     n_final_iters = 5)  # this is not enough
#' # plot clusters:
#' scols = c(cellcols, brewer.pal(12, "Set3")[-2], brewer.pal(8, "Set2"))[1:length(unique(semi$clust))]
#' names(scols) = unique(semi$clust)
#' plot(cellannot$x, cellannot$y, pch = 16, col = scols[semi$clust])
#' # draw flightpath plot for summarizing cell typing results:
#' fp = flightpath_layout(probs = semi$probs, cluster_xpos = NULL, cluster_ypos = NULL)
#' plot(fp$cellpos, pch = 16, col = alpha(scols[semi$clust], 0.5))
#' text(fp$clustpos[, 1], fp$clustpos[, 2], rownames(fp$clustpos), cex = 1.5)
cellEMClust <- function(counts, neg, bg = NULL, init_clust = NULL, n_clusts = NULL,
                        fixed_profiles = NULL, align_genes = TRUE, nb_size = 10, n_iters = 20,
                        method = "CEM", shrinkage = 0.8,
                        subset_size = 1000, n_starts = 10, n_benchmark_cells = 500,
                        n_final_iters = 3, pct_drop = 1/10000, min_prob_increase = 0.05) {

  # align genes in counts and fixed_profiles
  if (align_genes & !is.null(fixed_profiles)) {
    sharedgenes <- intersect(rownames(fixed_profiles), colnames(counts))
    lostgenes <- setdiff(colnames(counts), rownames(fixed_profiles))

    # subset:
    counts <- counts[, sharedgenes]
    fixed_profiles <- fixed_profiles[sharedgenes, ]
    if (is.matrix(bg)) {
      bg <- bg[, sharedgenes]
    }

    # warn about genes being lost:
    if ((length(lostgenes) > 0) & length(lostgenes < 50)) {
      message(paste0("The following genes in the count data are missing from fixed_profiles and will be omitted: ",
                     paste0(lostgenes, collapse = ",")))
    }
    if (length(lostgenes) > 50) {
      message(paste0(length(lostgenes), " genes in the count data are missing from fixed_profiles and will be omitted"))
    }

  }

  # infer bg if not provided: assume background is proportional to the scaling factor s
  if (is.null(bg)) {
    s <- Matrix::rowMeans(counts)
    bgmod <- stats::lm(neg ~ s - 1)
    bg <- bgmod$fitted
  }

  # define random starts:
  randstarts <- lapply(seq_len(n_starts), function(i){
    sample(x = seq_len(nrow(counts)), size = subset_size, replace = FALSE)
  })

  # select subset of cells for comparing the random starts:
  benchmarking_cells <- sample(x = seq_len(nrow(counts)), size = n_benchmark_cells, replace = FALSE)

  # loglik for comparing separate starts:
  best_loglik <- -Inf
  best_clust <- NULL
  # cluster from each random start:
  for (i in seq_len(n_starts)) {
    message(paste0("random start ", i))
    use = randstarts[[i]]

    if (is.null(init_clust)) {
      randinit <- NULL
    } else {
      randinit <- init_clust[use]
    }
    # run clustering:
    tempclust <- nbclust(counts = counts[use, ],
                         neg = neg[use], bg = bg[use],
                         init_clust = randinit, n_clusts = n_clusts,
                         fixed_profiles = fixed_profiles,
                         nb_size = nb_size, n_iters = n_iters,
                         method = method, shrinkage = shrinkage,
                         updated_reference = NULL,
                         pct_drop = pct_drop, 
                         min_prob_increase = min_prob_increase)

    # now get the loglik of the benchmarking cells under this clustering scheme:
    loglik_thisclust <- apply(tempclust$profiles, 2, function(ref) {
      lldist(x = ref,
             mat = counts[benchmarking_cells, ],
             bg = bg[benchmarking_cells],
             size = nb_size)
    })
    mean_loglik_this_clust <- mean(apply(loglik_thisclust, 1, max))
    # save this clustering if it's the best so far:
    if (mean_loglik_this_clust > best_loglik) {
      best_clust <- tempclust
      best_loglik <- mean_loglik_this_clust
    }
  } # now on to the final clustering

  if (is.null(init_clust)) {
    # for the final clustering, get initial cell type assignments using the best subset clustering result:
    logliks_under_init_clust <- apply(best_clust$profiles, 2, function(ref) {
      lldist(x = ref, mat = counts, bg = bg, size = nb_size)
    })
    final_clust_init <- colnames(logliks_under_init_clust)[
      apply(logliks_under_init_clust, 1, which.max)]
  }
  else {
    final_clust_init <- init_clust
  }

  # now run the final clustering:
  message("clustering all cells")
  finalclust <- nbclust(counts = counts, neg = neg, bg = bg,
                        init_clust = final_clust_init, n_clusts = 0,
                        fixed_profiles = fixed_profiles, nb_size = nb_size,
                        n_iters = n_final_iters,
                        method = method, shrinkage = shrinkage,
                        updated_reference = best_clust$updated_reference,
                        pct_drop = pct_drop,
                        min_prob_increase = min_prob_increase)

  return(finalclust)
}
