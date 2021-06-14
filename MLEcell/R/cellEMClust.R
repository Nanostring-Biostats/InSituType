
# function for calculating loglikelihood distance between cell profiles (mat) 
# and a reference matrix (x).
# (duplicate with the same function found in findCellTypes, so commenting out here)
lldist <- function(x, mat, bg = 0.01, size = 10, digits = 2) {
  
  # convert to matrix form if only a vector was input:
  if (is.vector(mat)) {
    mat = matrix(mat, nrow = 1)
  }
  
  # calc scaling factor to put y on the scale of x:
  bgsub = pmax(sweep(mat, 1, bg, "-"), 0)
  s = Matrix::rowSums(bgsub) / sum(x) 
  # override it if s is negative:
  s[s <= 0] = Matrix::rowSums(mat[s <= 0, , drop = FALSE]) / sum(x)
  
  # expected counts:
  yhat = sweep(s %*% t(x), 1, bg, "+")
  
  # loglik:
  lls = dnbinom(x = as.matrix(mat), size = size, mu = yhat, log = T)
  
  return(round(rowSums(lls), digits))
}


#' M step
#' 
#' Compute probability that each cell belongs to a given cluster
#' @param counts Counts matrix, cells * genes.
#' @param means Matrix of mean cluster profiles,
#'  with genes in rown and clusters in columns.
#' @param s Vector of scaling factors for each cell, 
#'   e.g. as defined by cell area. 
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
  logliks = sweep(logliks, 1, apply(logliks, 1, max), "-")
  # get on likelihood scale:
  liks = exp(logliks)
  # convert to probs
  probs = sweep(liks, 1, rowSums(liks), "/")
  
  return(round(probs, digits))
}


#' E step: estimate each cluster's mean profile
#'
#' Given cell assignments (or posterior probabilities), estimate the mean 
#'  profile of each cluster.
#' @param counts Counts matrix, cells * genes.
#' @param clust Vector of cluster assignments, or a matrix of probabilities 
#'   of cells (rows) belonging to clusters (columns).
#' @param s Vector of scaling factors for each cell, 
#'   e.g. as defined by cell area. 
#' @param neg Vector of mean background counts
#' @return A matrix of cluster profiles, genes * clusters
Estep <- function(counts, clust, s, neg) {
  

  # get cluster means:
  if (is.vector(clust)) {
    means = sapply(unique(clust), function(cl) {
      pmax(colSums(counts[clust == cl, , drop = FALSE]) - sum(neg[clust == cl]), 0)
    })
  }
  if (is.matrix(clust)) {
    means <- apply(clust, 2, function(x) {
      wts <- x / sum(x)
      pmax(colSums(sweep(counts, 1, wts, "*")) - sum(neg * wts), 0)
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
#' @param s Vector of scaling factors for each cell, 
#'   e.g. as defined by cell area. 
#' @param neg Vector of mean background counts
#' @param fixed_profiles Matrix of expression profiles of pre-defined clusters,
#'  e.g. from previous scRNA-seq. These profiles will not be updated by the EM algorithm.
#'  Colnames must all be included in the init_clust variable.
#' @param shrinkage Fraction by which to shrink the average profiles towards
#'  the fixed profiles. 1 = keep the fixed profile; 0 = don't shrink the mean profile.
#' @return A matrix of cluster profiles, genes * clusters
Estep_reference <- function(counts, clust, s, neg, fixed_profiles, shrinkage = 0.5) {
  

  # get cluster means:
  if (is.vector(clust)) {
    means = sapply(unique(clust), function(cl) {
      pmax(colSums(counts[clust == cl, , drop = FALSE]) - sum(neg[clust == cl]), 0)
    })
  }
  if (is.matrix(clust)) {
    means <- apply(clust, 2, function(x) {
      wts <- x / sum(x)
      pmax(colSums(sweep(counts, 1, wts, "*")) - sum(neg * wts), 0)
    })
  }
  # scale the fixed_profiles to be on the same scale as the means:
  scaled_fixed_profiles <- means * NA
  for (name in colnames(means)) {
    scaled_fixed_profiles[, name] = 
      fixed_profiles[, name] * mean(means[, name]) / mean(fixed_profiles[, name])
  }
  
  # shrink means towards the fixed profiles:
  updated_profiles <- means * (1 - shrinkage) + scaled_fixed_profiles * shrinkage
  
  return(updated_profiles)
}


#' Estimate NB size parameter
#' 
#' Given cell assignments (or posterior probabilities), update theta. 
#' Assumption is one theta for all data. 
#' @param counts Counts matrix, cells * genes.
#' @param clust Vector of cluster assignments, or a matrix of probabilities 
#'   of cells (rows) belonging to clusters (columns).
#' @param s Vector of scaling factors for each cell, 
#'   e.g. as defined by cell area. 
#' @param bg Expected background
#' @param means Matrix of means
#' @return A scalar giving the mle for the size parameter
Estep_size <- function(counts, clust, s, bg) {
  
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

#' Cluster via EM alogorithm based on cell logliks
#' 
#' Cluster single cell gene expression data using an EM algorithm.
#' @param counts Counts matrix, cells * genes.
#' @param s Vector of scaling factors for each cell, 
#'   e.g. as defined by cell area. 
#' @param neg Vector of mean negprobe counts per cell   
#' @param bg Expected background
#' @param init_clust Vector of initial cluster assignments. 
#' If NULL, initial assignments will be automatically inferred.
#' @param n_clusts Number of clusters, in addition to any pre-specified cell types. 
#' @param fixed_profiles Matrix of expression profiles of pre-defined clusters,
#'  e.g. from previous scRNA-seq. These profiles will not be updated by the EM algorithm.
#'  Colnames must all be included in the init_clust variable.
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
#' @return A list, with the following elements:
#' \enumerate{
#' \item probs: a matrix of probabilies of all cells (rows) belonging to all clusters (columns)
#' \item profiles: a matrix of cluster-specific expression profiles
#' }
nbclust <- function(counts, s, neg, bg = NULL, init_clust = NULL, n_clusts = NULL,
                    fixed_profiles = NULL, nb_size = 10, n_iters = 20, 
                    method = "CEM", shrinkage = 0.8, subset_size = 1000) {
  
  # checks:
  if (min(s) <= 0) {
    stop("scaling factors must be positive numbers.")
  }
  
  # infer bg if not provided: assume background is proportional to the scaling factor s
  if (is.null(bg)) {
    bgmod <- lm(neg ~ s - 1)
    bg <- bgmod$fitted
  }
  
  # specify which profiles to leave fixed:
  keep_profiles = NULL
  if (!is.null(fixed_profiles)) {
    keep_profiles = colnames(fixed_profiles)
  }
  
  ### get initial profiles:
  
  # if no initial clustering is available, quickly learn profiles:
  if (is.null(init_clust)) {
    if (is.null(n_clusts)) {
      stop("Must specify either init_clust or n_clusts.")
    }
    n_fixed_profiles = 0
    if (!is.null(fixed_profiles)) {
      n_fixed_profiles = ncol(fixed_profiles)
    }
    clustnames = c(colnames(fixed_profiles), letters, paste0(rep(letters, each = 26), letters))[
      seq_len(n_clusts + n_fixed_profiles)]
    # arbitrary but non-random initialization:
    init_clust = rep(clustnames, ceiling(nrow(counts) / length(clustnames)))[seq_len(nrow(counts))]
  }
  
  
  # for deriving initial profiles, subset on only the cells that aren't part of a pre-specified cluster:
  tempuse = !is.element(init_clust, keep_profiles)
  bgtemp = bg
  if (length(bg) == nrow(counts)) {
    bgtemp = bg[tempuse]
  }
  if (is.matrix(bg)) {
    bgtemp = bg[tempuse, ]
  }
  
  # if an initial clustering is available, use it to estimate initial profiles:
  if (!is.null(init_clust)) {
    new_profiles <- Estep(counts = counts[tempuse, ], 
                          clust = init_clust[tempuse], 
                          s = s[tempuse], neg = neg[tempuse]) 
    profiles <- cbind(fixed_profiles, new_profiles)
  }
  
  clust_old = init_clust
  n_changed = c()
  
  for (iter in seq_len(n_iters)) {
    message(paste0("iter ", iter))
    # M-step: get cell * cluster probs:
    probs <- Mstep(counts = counts, 
                   means = profiles, 
                   bg = bg, 
                   size = nb_size) 
    
    # E-step: update profiles:
    if (method == "CEM") {
      tempprobs = probs[, setdiff(colnames(probs), keep_profiles)]
      # update the new cluster profiles:
      new_profiles <- Estep(counts = counts, 
                            clust = tempprobs,
                            s = s, 
                            neg = neg)
      # update the reference profiles / "fixed_profiles" \
      updated_reference <- 
        Estep_reference(counts = counts, 
                        clust = probs[, colnames(fixed_profiles)],
                        s = s, 
                        neg = neg,
                        fixed_profiles = fixed_profiles,
                        shrinkage = shrinkage)
      #nb_size <- Estep_size(counts = counts, )
      
    }
    if (method == "EM") {
      tempprobs = 1 * t(apply(probs[, setdiff(colnames(probs), keep_profiles)], 1, ismax))
      # update the new cluster profiles:  
      new_profiles <- Estep(counts = counts, 
                            clust = tempprobs,
                            s = s, 
                            neg = neg)
      # update the reference profiles / "fixed_profiles" 
      updated_reference <- 
        Estep_reference(counts = counts, 
                        clust = probs[, colnames(fixed_profiles)],
                        s = s, 
                        neg = neg,
                        fixed_profiles = fixed_profiles,
                        shrinkage = shrinkage)
      
      # for any profiles that have been lost, replace them with their previous version:
      lostprofiles = names(which(colSums(!is.na(new_profiles)) == 0))
      new_profiles[, lostprofiles] = profiles[, lostprofiles]
    }
    
    
    profiles <- cbind(updated_reference, new_profiles)

    # get cluster assignment
    clust = colnames(probs)[apply(probs, 1, which.max)]
    
    # record number of switches
    n_changed = c(n_changed, sum(clust != clust_old))
    clust_old = clust
  }
  
  # get loglik of each cell:
  logliks = unlist(sapply(seq_len(nrow(counts)), function(i) {
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
             n_changed = n_changed,
             logliks = logliks)
  return(out)
}

#' For a numeric object, return a logical object of whether each element is the max or not. 
ismax <- function(x) {
  return(x == max(x, na.rm = T))
}



#' Clustering wrapper function
#' 
#' A wrapper for nbclust, to manage subsampling and multiple random starts
#' @param counts Counts matrix, cells * genes.
#' @param s Vector of scaling factors for each cell, 
#'   e.g. as defined by cell area. 
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
#' @return A list, with the following elements:
#' \enumerate{
#' \item probs: a matrix of probabilies of all cells (rows) belonging to all clusters (columns)
#' \item profiles: a matrix of cluster-specific expression profiles
#' }
#' @export
#' @examples
#' # load data ("raw" and "cellannot"):
#' load("inst/extdata/concise melanoma 300plex data.RData") 
#' # predict per-cell bbackground:
#' bgmodel = lm(rowSums(raw[, grepl("NegPrb", colnames(raw))]) ~ rowSums(raw) - 1)$coef
#' bg.predicted = rowSums(raw) * bgmodel
# run unsupervised clustering with several random starts:
#' unsup <- cellEMClust(counts = raw, 
#'                      s = pmax(rowSums(raw), 1), 
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
#'                     s = pmax(rowSums(raw), 1), 
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
cellEMClust <- function(counts, s, neg, bg = NULL, init_clust = NULL, n_clusts = NULL,
                        fixed_profiles = NULL, align_genes = TRUE, nb_size = 10, n_iters = 20, 
                        method = "CEM", shrinkage = 0.8, 
                        subset_size = 1000, n_starts = 10, n_benchmark_cells = 500,
                        n_final_iters = 3) {
  
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
    bgmod <- lm(neg ~ s - 1)
    bg <- bgmod$fitted
  }
  
  # define random starts:
  randstarts = vector("list", n_starts)
  for (i in seq_len(n_starts)) {
    randstarts[[i]] <- sample(x = seq_len(nrow(counts)), size = subset_size, replace = FALSE)
  }
  
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
    tempclust <- nbclust(counts = counts[use, ], s = s[use], 
                         neg = neg[use], bg = bg[use], 
                         init_clust = randinit, n_clusts = n_clusts,
                         fixed_profiles = fixed_profiles, 
                         nb_size = nb_size, n_iters = n_iters, 
                         method = method, shrinkage = shrinkage)
    
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
    logliks_under_init_clust <- apply(tempclust$profiles, 2, function(ref) {
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
  finalclust <- nbclust(counts = counts, s = s, neg = neg, bg = bg, 
                        init_clust = final_clust_init, n_clusts = NULL,
                        fixed_profiles = fixed_profiles, nb_size = nb_size, 
                        n_iters = n_final_iters, 
                        method = method, shrinkage = shrinkage)
  
  return(finalclust)
}