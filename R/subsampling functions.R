# functions for handling the subsampling operations of the main clustering function



#' Clustering function with 4 levels of subsampling
#' 
#' @param sketchingdata Optional matrix of data for use in non-random sampling via "sketching".
#'  If not provided, then the data's first 20 PCs will be used. 
insitutype <- function(counts, neg, bg = NULL, 
                       n_clusts = NULL,
                       fixed_profiles = NULL, 
                       sketchingdata = NULL,
                       align_genes = TRUE, nb_size = 10, 
                       method = "CEM", 
                       init_clust = NULL, n_starts = 10, n_benchmark_cells = 50000,
                       n_phase1 = 5000, n_phase2 = 20000, n_phase3 = 100000,
                       pct_drop = 1/10000, min_prob_increase = 0.05)
  
  #### preliminaries ---------------------------
  
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

  # get data for subsetting if not already provided
  # (e.g., if PCA is the choice, then point to existing PCA results, and run PCA if not available
  if (!is.null(sketchingdata) {
    # check that it's correct:
    if (nrow(sketchingdata) != nrow(counts)) {
      warning("counts and sketchingdata have different numbers of row. Discarding sketchingdata.")
      sketchingdata <- NULL
    }
  }
  if (is.null(sketchingdata)) {
    sketchingdata <- prepDataForSketching(counts)
  }
  n_phase1 = min(n_phase1, nrow(counts))
  n_phase2 = min(n_phase2, nrow(counts))
  n_phase3 = min(n_phase3, nrow(counts))
  n_benchmark_cells = min(n_benchmark_cells, nrow(counts))
  
  
  #### phase 1: many random starts in small subsets -----------------------------
  
  if (!is.null(init_clust)) {
    message("init_clust was provided, so phase 1 - random starts in small subsets - will be skipped.")
    
    tempprofiles <- sapply(by(counts[!is.na(init_clust), ], init_clust[!is.na(init_clust)], colMeans), cbind)
    rownames(tempprofiles) <- colnames(counts)

  }
  else {
    message(paste0("phase 1: random starts in ", n_phase1, " cell subsets"))
    
    # get a list in which each element is the cell IDs to be used in a subset
    random_start_subsets <- list()
    for (i in 1:n_starts) {
      random_start_subsets[[i]] <- geoSketch(X = sketchingdata,
                                            N = n_phase1,
                                            alpha=0.1,
                                            max_iter=200,
                                            returnBins=FALSE,
                                            minCellsPerBin = 1,
                                            seed=NULL)#$sampledCells
      # NOTE: should probably run plaid calculation just once, then sample across plaids multiple times
      #random_start_subsets[[i]] <- sample(rownames(counts), n_phase1, replace = F)
    }
    
    # get a vector of cells IDs to be used in comparing the random starts:
    #benchmarking_subset <- sample(rownames(counts), n_benchmark_cells, replace = F)
    benchmarking_subset <- geoSketch(X = sketchingdata,
                             N = n_benchmark_cells,
                             alpha=0.1,
                             max_iter=200,
                             returnBins=FALSE,
                             minCellsPerBin = 1,
                             seed=NULL)
    
    # run nbclust from each of the random subsets, and save the profiles:
    profiles_from_random_starts <- list()
    for (i in 1:n_starts) {
      profiles_from_random_starts[[i]] <- nbclust(
        counts = counts[random_start_subsets[[i]], ], 
        neg = neg[random_start_subsets[[i]]], 
        bg = bg[random_start_subsets[[i]]],
        init_clust = NULL, 
        n_clusts = n_clusts,
        fixed_profiles = fixed_profiles, 
        nb_size = nb_size,
        method = method, 
        updated_reference = NULL,
        pct_drop = pct_drop,
        min_prob_increase = min_prob_increase
      )$profiles
    }
    
    # find which profile matrix does best in the benchmarking subset:
    benchmarking_logliks <- c()
    for (i in 1:n_starts) {
      templogliks <- apply(profiles_from_random_starts[[i]], 2, function(ref) {
        lldist(x = ref,
               mat = counts[benchmarking_subset, ],
               bg = bg[benchmarking_subset],
               size = nb_size)
      })
      # take the sum of cells' best logliks:
      benchmarking_logliks[i] = sum(apply(templogliks, 1, max))
    }
    best_start <- which.max(benchmarking_logliks)
    tempprofiles <- profiles_from_random_starts[[best_start]]
    
    rm(profiles_from_random_starts)
    rm(templogliks)
  }
  
  #### phase 2: -----------------------------------------------------------------
  message(paste0("phase 2: refining best random start in a ", n_phase2, " cell subset"))
  phase2_sample <- geoSketch(X = sketchingdata,
                             N = n_phase2,
                             alpha=0.1,
                             max_iter=200,
                             returnBins=FALSE,
                             minCellsPerBin = 1,
                             seed=NULL)
  
  # get initial cell type assignments:
  if (!is.null(init_clust)) {
    temp_init_clust = init_clust[phase2_sample]
  } else {
    templogliks <- apply(tempprofiles, 2, function(ref) {
      lldist(x = ref,
             mat = counts[phase2_sample, ],
             bg = bg[phase2_sample],
             size = nb_size)
    })
    temp_init_clust <- colnames(templogliks)[apply(templogliks, 1, which.max)]
    rm(templogliks)
  }
  
  # run nbclust, initialized with the cell type assignments derived from the previous phase's profiles
  clust2 <- nbclust(counts = counts[phase2_sample, ], 
                    neg = neg[phase2_sample], 
                    bg = bg[phase2_sample],
                    init_clust = temp_init_clust, 
                    n_clusts = n_clusts,
                    fixed_profiles = fixed_profiles, 
                    nb_size = nb_size,
                    method = method, 
                    updated_reference = NULL,   #<---- for now, decision is to not use updated_reference from phase1 due to instability arising from small n.
                    pct_drop = pct_drop,
                    min_prob_increase = min_prob_increase)
  
  #<< what should we save from this step? >>
  tempprofiles <- clust2$profiles
  #if (!is.null(fixed_profiles)) {
  #  tempprofiles[, colnames(fixed_profiles), drop = FALSE] <- fixed_profiles[rownames(tempprofiles), ]
  #}   # <----- commented out this part since tempprofiles are used to get init_clusts in next round, so should benefit from platform effects estimation
  
  #### phase 3: -----------------------------------------------------------------
  message(paste0("phase 3: finalizing clusters in a ", n_phase3, " cell subset"))
  
  #phase3_sample <- geoSketch(X = get(sketchingdataname),
  #                           N = n_phase3,
  #                           alpha=0.1,
  #                           max_iter=200,
  #                           returnBins=FALSE,
  #                           minCellsPerBin = 1,
  #                           seed=NULL)$sampledCells
  phase3_sample <- sample(rownames(counts), n_phase3)
  
  # get initial cell type assignments:
  templogliks <- apply(tempprofiles, 2, function(ref) {
    lldist(x = ref,
           mat = counts[phase3_sample, ],
           bg = bg[phase3_sample],
           size = nb_size)
  })
  temp_init_clust <- colnames(templogliks)[apply(templogliks, 1, which.max)]
  rm(templogliks)
  
  # run nbclust, initialized with the cell type assignments derived from the previous phase's profiles
  clust3 <- nbclust(counts = counts[phase3_sample, ], 
                    neg = neg[phase3_sample], 
                    bg = bg[phase3_sample],
                    init_clust = temp_init_clust, 
                    n_clusts = 0,
                    fixed_profiles = fixed_profiles, 
                    nb_size = nb_size,
                    method = method, 
                    updated_reference = clust2$updated_reference, 
                    pct_drop = pct_drop,
                    min_prob_increase = min_prob_increase)
  
  #<< what should we save from this step? >>
  # (copy from previous)
  
  #### phase 4: -----------------------------------------------------------------
  message(paste0("phase 4: classifying all ", nrow(counts), " cells"))
  
  profiles <- clust3$profiles
  
  logliks <- apply(profiles, 2, function(ref) {
    lldist(x = ref,
           mat = counts,
           bg = bg,
           size = nb_size)
  })
  clust <- colnames(logliks)[apply(logliks, 1, which.max)]
  
  templogliks <- sweep(logliks, 1, apply(logliks, 1, max ), "-" )
  # get on likelihood scale:
  liks <- exp(templogliks)
  # convert to probs
  probs <- sweep(liks, 1, rowSums(liks), "/")
  
  
  out = list(clust = clust,
             probs = round(probs, 3),
             profiles = sweep(profiles, 2, colSums(profiles), "/") * nrow(profiles),
             logliks = round(logliks, 3))
  
}



#' Prepare data for geoSketch 
#' 
#' Process raw counts data for input into geoSketching. 
#' @param counts Counts matrix
#' @param method What kind of data to extract. 
#' @return A matrix of data for geoSketch, with cells in rows and features in columns
#' @importFrom irlba prcomp_irlba
prepDataForSketching <- function(counts) {
  # get PCs:
  scaling_factors <- pmax(apply(counts, 2, quantile, 0.99), 5)
  pcres <- irlba::prcomp_irlba(x = sweep(counts, 2, scaling_factors), n = 20, retx = TRUE, center = TRUE, scale. = FALSE)$x
  rownames(pcres) <- rownames(counts)
  return(pcres)
}






#' Function for creating a biased sample of a given dataset with the aim of retaining cells with unique expression vectors
#'
#' @export geoSketch


#' @examples
#' 
#' # Example use of geometric sketching to return sampled list of cellIDs and binIDs
#' library(Ptolemy) # Load library
#' data(mini_tma) # Load data
#' 
#' mini_tma <- runPCA(mini_tma, expression_values="raw", ncp = 5) # Run PCA
#' X <- Giotto:::select_dimReduction(mini_tma, name = "pca") # Use PCA results as expression matrix
#' N <- 100 # Select the desired sample size
#' 
#' sampledCells <- geoSketch(X, N) # Generate list of sampled cells
#' geoSketchRes <- geoSketch(X, N, returnBins=TRUE) # Generate list of sampled cells and return named vector of binIDs used for sampling


#' @title geoSketch
#' Function for sampling cells evenly across expression space
#'
#' Given an expression matrix, evenly bin cells across expression space and return bin labels and/or sampled cellIDs
#' @param X feature matrix with cellIDs as rows and featureIDs as columns (can be counts, normalized expression, PCA, UMAP, etc.)
#' @param N desired sample size
#' @param alpha defines the acceptable minimum number of bins to sample from as `(1-alpha)*N`
#' @param max_iter maximum number of iterations used to achieve an acceptable minimum number of bins
#' @param returnBins determines whether or not to pass back bin labels for each cell
#' @param minCellsPerBin the minimum number of cells required for a bin to be considered for sampling
#' @param seed set seed for random sampling
#' 
#' @return sampledCells, a vector of cellIDs sampled using the geometric sketching method
#' @return Plaid, a named vector of binIDs where names correspond to cellIDs
geoSketch <- function(X, N,
                      alpha=0.1,
                      max_iter=200,
                      returnBins=FALSE,
                      minCellsPerBin = 1,
                      seed=NULL){
  
  # Define seed for sampling if given
  if (!is.null(seed)){
    set.seed(seed)
  }
  
  # Determine the total number of cells and compare it to the desired sample size 
  nCells <- length(rownames(X))
  
  if (N > nCells){
    # Stop the function and return and error is the desired sample size is greater than the total number of cells given
    stop(paste0("N = ",N ," is greater than the number of cells in the feature matrix, nCells = ", nCells))
  }
  
  # Normalize features on a range between 0 and 1 in order to facilitate even binning across expression space
  X <- apply(X, 2, function(Y) (Y-min(Y))/(max(Y)-min(Y)))
  
  # Iterate the number of bins each feature is broken into until the total number of bins containing cells is larger than `(1-alpha)*N`
  iter <- 1 # Set starting iteration
  while (iter <= max_iter){ # Break loop when the number of iterations surpasses the defined maximum
    message(paste0("Iteration Number ", as.character(iter))) # Report the current iteration the loop is on
    
    bins <- seq(0, 1, length=2+iter) # Define the bin ranges based on the current iteration
    Xbins <- apply(X, 2, function(Y) .bincode(Y, bins, TRUE, TRUE)) # Bin cells across each feature using the given bin ranges
    
    Plaid <- apply(Xbins, 1, paste, collapse="") # Collapse the bin assignments to form a multidimensional binID
    
    PlaidCounts <- table(Plaid) # Count the number of cells per bin
    PlaidMinCells <- names(PlaidCounts[PlaidCounts >= minCellsPerBin]) # Grab binIDs for bins that contain the minimum number of cells or more
    
    nPlaid <- length(PlaidMinCells) # Determine the number of unique bins containing at least the minimum number of cells
    
    # Check to see if the total number of bins containing cells is larger than `(1-alpha)*N`
    if (nPlaid<(1-alpha)*N){
      iter <- iter + 1 # If not continue onto the the next iteration
    } else {
      break # If the desired minimum number of bins has been reached end the loop
    }
    
    nCellsInBins <- length(Plaid[Plaid %in% PlaidMinCells])
    if (N > nCellsInBins){
      # Stop the function and return and error is the desired sample size is greater than the total number of cells in the bins to be sampled from
      stop(paste0("N = ",N ," is greater than the number of cells in bins above the minimum size, nCellsInBins = ", nCellsInBins,
                  "\n  Please modify the parameters N, minCellsPerBin, and/or alpha before attempting to resample"))
    }
    
  }
  
  # Determine and report back the average number of cells assigned to each unique bin
  # This is a useful metric for understanding how distinct the geometric sketch will be from random sampling
  # If the number of cells assigned to each unique bin is 1 then the geometric sketch is equivalent to random sampling
  cellsPerPlaid <- round(mean(table(Plaid)),2)
  message(paste0(cellsPerPlaid, " cells per geometric bin."))
  
  sampledCells = c() # Create the vector which sampled cells will be added to
  names(Plaid) <- row.names(X) # Ensure binIDs are associated with their respective cells
  
  # Sample only from the cells left in bins with more than the the minimum number of cells
  PlaidLeftover <- Plaid[Plaid %in% PlaidMinCells]
  
  sampledCells <- c() # Begin sampling cells
  while (length(sampledCells) < N){
    
    # After each round of sampling ensure additional sampling is performed on unsampled cells only
    PlaidLeftover <- PlaidLeftover[!(names(PlaidLeftover) %in% sampledCells)]
    PlaidAddrsLeftover <- unique(PlaidLeftover)
    
    nCellsLeft <- N - (length(sampledCells)) # Calculate the number of cells left to sample
    if (nCellsLeft > length(PlaidAddrsLeftover)){ # If the number of cells left to sample is greater than the number of bins left, sample one cell from each bin
      sampleReamaining <- sapply(PlaidAddrsLeftover,
                                 FUN = function(PlaidAddr){names(sample(PlaidLeftover[PlaidLeftover == PlaidAddr], 1))})
    } else{ # Else sample the required number of cells one at a time from random remaining bins
      sampleReamaining <- sapply(sample(PlaidAddrsLeftover)[1:nCellsLeft],
                                 FUN = function(PlaidAddr){names(sample(PlaidLeftover[PlaidLeftover == PlaidAddr], 1))})
    }
    sampledCells <- c(sampledCells, sampleReamaining) # Add the newly sampled cells to the existing list of sampled cells
  }
  
  if (length(sampledCells) > N){ # In the event where the number of Bins was greater than N, randomly sample N cells from the already sampled cells
    sampledCells <- sample(sampledCells, N)
  }
  
  # Return sampled cellIDs and also cell binIDs if requested
  if (!returnBins){
    return(sampledCells)
  } else{
    return(list(sampledCells = sampledCells, binIDs = Plaid))
  }
}
