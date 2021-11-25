
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
#' @param freq a vector of cells frequencies summing up equal to 1. 
#' @param bg Expected background
#' @param size NB size parameter
#' @param digits Round the output to this many digits (saves memory)
#' @param return_loglik If TRUE, logliks will be returned. If FALSE, probabilities will be returned. 
#' @return Matrix of probabilities of each cell belonging to each cluster
#' @export
Mstep <- function(counts, means, freq, bg = 0.01, size = 10, digits = 2, return_loglik = FALSE) {
  # get logliks of cells * clusters
  logliks <- apply(means, 2, function(x) {
    lldist(x = x, mat = counts, bg = bg, size = size)
  })
  if (return_loglik) {
    return(round(logliks, digits))
  } else {
    # first rescale (ie recenter on log scale) to avoid rounding errors:
    logliks <- sweep( logliks , 1 , apply( logliks , 1 , max ) , "-" )
    # get on likelihood scale:
    liks <- sweep(exp( logliks ), 2, freq, "*")
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
#'
#' @export
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




#' Estimate NB size parameter
#'
#' Given cell assignments (or posterior probabilities), update theta.
#' Assumption is one theta for all data.
#' @param counts Counts matrix, cells * genes.
#' @param clust Vector of cluster assignments, or a matrix of probabilities
#'   of cells (rows) belonging to clusters (columns).
#' @param bg Expected background
#' @return A scalar giving the mle for the size parameter
#' @export
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
#' @param anchors Vector giving "anchor" cell types, for use in semi-supervised clustering. 
#'  Vector elements will be mainly NA's (for non-anchored cells) and cell type names
#'  for cells to be held constant throughout iterations. 
#' @param init_profiles Matrix of cluster profiles under which to begin iterations.
#' If NULL, initial assignments will be automatically inferred, using init_clust 
#' if available, and using random clusters if not. 
#' @param init_clust Vector of initial cluster assignments.
#' If NULL, initial assignments will be automatically inferred.
#' @param n_clusts Number of clusters, in addition to any pre-specified cell types.
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param  Numer of iterations
#' @param method Whether to run a SEM algorithm (points given a probability
#'   of being in each cluster) or a classic EM algorithm (all points wholly
#'   assigned to one cluster).
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
nbclust <- function(counts, neg, bg = NULL, anchors = NULL,
                    init_profiles = NULL, init_clust = NULL, n_clusts = NULL,
                    nb_size = 10,
                    method = "CEM", 
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

  clusterlog = NULL

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
    profiles <- Estep(counts = counts,
                           clust = init_clust,
                           neg = neg)
    clustnames <- unique(init_clust)
  }
  
  #### run EM algorithm iterations: ----------------------------------
  pct_changed = c()
  if (logresults) {
    clusterlog <- init_clust
  }
  
  profiles_freq <- setNames(rep(1/ncol(profiles), ncol(profiles)), colnames(profiles))

  for (iter in seq_len(max_iters)) {
    message(paste0("iter ", iter))
    # M-step: get cell * cluster probs:
    
    probs <- Mstep(counts = counts,
                   means = profiles,
                   freq = profiles_freq, 
                   bg = bg,
                   size = nb_size)
    # override assignments for anchor cells
    if (!is.null(anchors)) {
      for (cell in setdiff(unique(anchors), NA)) {
        temprows <- which((anchors == cell) & !is.na(anchors))
        probs[temprows, cell] <- 1 #rep(1, sum((anchors == cell) & !is.na(anchors)))
        probs[temprows, setdiff(colnames(probs), cell)] <- 0
      }
    }
    if (logresults) {
      clusterlog <- cbind(clusterlog, colnames(probs)[apply(probs, 1, which.max)])
    }
    
    oldprofiles <- profiles

    # E-step: update profiles:
    if (method == "CEM") {
      tempprobs = probs
      # update the new cluster profiles:
      if (n_clusts != 0){
        profiles <- Estep(counts = counts,
                          clust = tempprobs,
                          neg = neg)
      }
    }
    
    if (method == "EM") {
      # update the new cluster profiles:
      if (n_clusts != 0){
        tempprobs = 1 * t(apply(probs, 1, ismax))
        profiles <- Estep(counts = counts,
                          clust = tempprobs,
                          neg = neg)        
      }
    }
    
    # for any profiles that have been lost, replace them with their previous version:
    lostprofiles = names(which(colSums(!is.na(profiles)) == 0))
    profiles[, lostprofiles] = oldprofiles[, lostprofiles]
    
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
    profiles_freq <- setNames(data.frame(prop.table(table(clust_old)))$Freq, 
                              data.frame(prop.table(table(clust_old)))$clust_old)
    profiles_freq <- profiles_freq[colnames(profiles)]
    # prevent any freqs from going all the way to zero:
    profiles_freq <- pmax(profiles_freq, 1e-3)
    probs_old_max = apply(probs, 1, max)
  }
  names(pct_changed) <- paste0("Iter_", seq_len(iter))
  ## get loglik of each cell:
  #logliks = unlist(sapply(seq_len(nrow(counts)), function(i) {   #<-------------------- should just assume vector background
  #  #if (length(bg) == 1) {
  #  #  bgtemp = bg
  #  #}
  #  #if (length(bg) == nrow(counts)) {
  #  #  bgtemp = bg[i]
  #  #}
  #  #if (is.matrix(bg)) {
  #  #  bgtemp = bg[i, ]
  #  #}
  #  bgtemp = bg[i]  # <--- now assuming vector-form background
  #  lldist(mat = counts[i, ], x = profiles[, clust[i]], bg = bgtemp, size = nb_size)
  #}))


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
#' @export
ismax <- function(x) {
  return(x == max(x, na.rm = T))
}


#' makeClusterNames
#' @param cNames existing cluster names
#' @param nClust number of clusters
#' @return a vector of complete cluster names
#'
#' @export
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

