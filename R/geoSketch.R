# functions for handling the subsampling operations of the main clustering function







#' Prepare data for geoSketch 
#' 
#' Process raw counts data for input into geoSketching. 
#' @param counts Counts matrix: cells x genes 
#' @param assay_type Assay type of RNA, protein (default = "rna")
#' 
#' @return A matrix of data for geoSketch, with cells in rows and features in columns
#' @importFrom irlba prcomp_irlba
#' @examples
#' data("mini_nsclc")
#' prepDataForSketching(counts=mini_nsclc$counts, assay_type="RNA")
prepDataForSketching <- function(counts, 
                                 assay_type = c("rna", "protein")) {
  assay_type <- match.arg(tolower(assay_type), c("rna", "protein"))
  
  # get PCs:
  if(identical(tolower(assay_type), "rna")){
    scaling_factors <- pmax(sparseMatrixStats::colQuantiles(counts, probs = 0.99), 5)
    x <- Matrix::t(Matrix::t(counts) / scaling_factors)
    pcres <- irlba::prcomp_irlba(x = x, n = min(20, ncol(counts) - 5), retx = TRUE, center = TRUE, scale. = FALSE)$x
  }
  if(identical(tolower(assay_type), "protein")){
    
    ## when the data is protein data 
    pcres <- irlba::prcomp_irlba(x = counts, n = min(20, ncol(counts) - 5), retx = TRUE, center = TRUE, scale. = TRUE)$x
    
  }
  rownames(pcres) <- rownames(counts)
  return(pcres)
}



#' Bin cells into "plaids"
#'
#' Assign cells to "plaids", very rough clusters.
#' @param X feature matrix with cellIDs as rows and featureIDs as columns (can
#'   be counts, normalized expression, PCA, UMAP, etc.)
#' @param N desired sample size
#' @param alpha defines the acceptable minimum number of bins to sample from as
#'   `(1-alpha)*N`
#' @param max_iter maximum number of iterations used to achieve an acceptable
#'   minimum number of bins
#' @param returnBins determines whether or not to pass back bin labels for each
#'   cell
#' @param minCellsPerBin the minimum number of cells required for a bin to be
#'   considered for sampling
#'
#' @return Plaid, a named vector of binIDs where names correspond to cellIDs
#' @examples
#' data("mini_nsclc")
#' geoSketch_get_plaid(mini_nsclc$counts, 100)
geoSketch_get_plaid <- function(X, N,
                                alpha=0.1,
                                max_iter=200,
                                returnBins=FALSE,
                                minCellsPerBin = 1) {
  
  # Determine the total number of cells and compare it to the desired sample size 
  nCells <- nrow(X)
  
  if (N > nCells) {
    # Stop the function and return and error is the desired sample size is greater than the total number of cells given
    stop(paste0("N = ", N, " is greater than the number of cells in the feature matrix, nCells = ", nCells))
  }
  
  # Normalize features on a range between 0 and 1 in order to facilitate even binning across expression space
  X <- apply(X, 2, function(Y) (Y - min(Y)) / (max(Y) - min(Y)))
  
  # Iterate the number of bins each feature is broken into until the total
  # number of bins containing cells is larger than `(1-alpha)*N`
  iter <- 1 # Set starting iteration
  while (iter <= max_iter) { # Break loop when the number of iterations surpasses the defined maximum
    #message(paste0("Iteration Number ", as.character(iter))) # Report the current iteration the loop is on
    
    bins <- seq(0, 1, length = 2 + iter) # Define the bin ranges based on the current iteration
    Xbins <- apply(X, 2, function(Y) .bincode(Y, bins, TRUE, TRUE)) # Bin cells across each feature using the given bin ranges
    
    Plaid <- apply(Xbins, 1, paste, collapse = "") # Collapse the bin assignments to form a multidimensional binID
    
    PlaidCounts <- table(Plaid) # Count the number of cells per bin
    PlaidMinCells <- names(PlaidCounts[PlaidCounts >= minCellsPerBin]) # Grab binIDs for bins that contain the minimum number of cells or more
    
    nPlaid <- length(PlaidMinCells) # Determine the number of unique bins containing at least the minimum number of cells
    
    # Check to see if the total number of bins containing cells is larger than `(1-alpha)*N`
    if (nPlaid < (1 - alpha) * N) {
      iter <- iter + 1 # If not continue onto the the next iteration
    } else {
      break # If the desired minimum number of bins has been reached end the loop
    }
    
    nCellsInBins <- length(Plaid[Plaid %in% PlaidMinCells])
    if (N > nCellsInBins) {
      # Stop the function and return and error is the desired sample size is greater than the total number of cells in the bins to be sampled from
      stop(
        paste0(
          "N = ",
          N,
          " is greater than the number of cells in bins above the minimum size, nCellsInBins = ",
          nCellsInBins,
          "\n  Please modify the parameters N, minCellsPerBin, and/or alpha before attempting to resample"
        )
      )
    }
    
  }
  
  # Determine and report back the average number of cells assigned to each unique bin
  # This is a useful metric for understanding how distinct the geometric sketch will be from random sampling
  # If the number of cells assigned to each unique bin is 1 then the geometric sketch is equivalent to random sampling
  cellsPerPlaid <- round(mean(table(Plaid)),2)
  message(paste0(cellsPerPlaid, " cells per geometric bin."))
  
  names(Plaid) <- row.names(X) # Ensure binIDs are associated with their respective cells
  
  return(Plaid)
}




#' Subsample from plaids
#' 
#' Sample cells, trying to give each plaid equal representation
#' @param Plaid Vector of cells' plaid IDs
#' @param N desired sample size
#' 
#' @return Plaid, a named vector of binIDs where names correspond to cellIDs
#' @return sampledCells, a vector of cellIDs sampled using the geometric sketching method
#' @examples
#' data("mini_nsclc")
#' plaids <- geoSketch_get_plaid(mini_nsclc$counts, 100)
#' geoSketch_sample_from_plaids(plaids, 5)
geoSketch_sample_from_plaids <- function(Plaid, N) {

  # define cells' sampling probabilities as the inverse of their plaid size:
  PlaidCounts <- table(Plaid) # Count the number of cells per bin
  prob <- 1 / PlaidCounts[Plaid]
  
  # subsample:
  sampledCells <- sample(seq_len(length(Plaid)), size = N, replace = FALSE, prob = prob)
  
  return(sampledCells)
}




#' Function for creating a biased sample of a given dataset with the aim of retaining cells with unique expression vectors
#'
#' @param X feature matrix with cellIDs as rows and featureIDs as columns (can be counts, normalized expression, PCA, UMAP, etc.)
#' @param N desired sample size
#' @param alpha defines the acceptable minimum number of bins to sample from as `(1-alpha)*N`
#' @param max_iter maximum number of iterations used to achieve an acceptable minimum number of bins
#' @param returnBins determines whether or not to pass back bin labels for each cell
#' @param minCellsPerBin the minimum number of cells required for a bin to be considered for sampling
#' 
#' @return sampledCells, a vector of cellIDs sampled using the geometric sketching method
#' @return Plaid, a named vector of binIDs where names correspond to cellIDs
#' @examples
#' data("mini_nsclc")
#' geoSketch(mini_nsclc$counts, 200)
geoSketch <- function(X, N,
                      alpha=0.1,
                      max_iter=200,
                      returnBins=FALSE,
                      minCellsPerBin = 1) {

  # Determine the total number of cells and compare it to the desired sample size 
  nCells <- nrow(X)
  
  if (N > nCells) {
    # Stop the function and return and error is the desired sample size is greater than the total number of cells given
    stop(paste0("N = ", N, " is greater than the number of cells in the feature matrix, nCells = ", nCells))
  }
  
  # Normalize features on a range between 0 and 1 in order to facilitate even binning across expression space
  X <- apply(X, 2, function(Y) (Y - min(Y)) / (max(Y) - min(Y)))
  
  # Iterate the number of bins each feature is broken into until the total number of bins containing cells is larger than `(1-alpha)*N`
  iter <- 1 # Set starting iteration
  while (iter <= max_iter) { # Break loop when the number of iterations surpasses the defined maximum
    bins <- seq(0, 1, length=2+iter) # Define the bin ranges based on the current iteration
    Xbins <- apply(X, 2, function(Y) .bincode(Y, bins, TRUE, TRUE)) # Bin cells across each feature using the given bin ranges
    
    Plaid <- apply(Xbins, 1, paste, collapse="") # Collapse the bin assignments to form a multidimensional binID
    
    PlaidCounts <- table(Plaid) # Count the number of cells per bin
    PlaidMinCells <- names(PlaidCounts[PlaidCounts >= minCellsPerBin]) # Grab binIDs for bins that contain the minimum number of cells or more
    
    nPlaid <- length(PlaidMinCells) # Determine the number of unique bins containing at least the minimum number of cells
    
    # Check to see if the total number of bins containing cells is larger than `(1-alpha)*N`
    if (nPlaid < (1 - alpha) * N) {
      iter <- iter + 1 # If not continue onto the the next iteration
    } else {
      break # If the desired minimum number of bins has been reached end the loop
    }
    
    nCellsInBins <- length(Plaid[Plaid %in% PlaidMinCells])
    if (N > nCellsInBins) {
      # Stop the function and return and error is the desired sample size is greater than the total number of cells in the bins to be sampled from
      stop(paste0("N = ", N, " is greater than the number of cells in bins above the minimum size, nCellsInBins = ", nCellsInBins,
                  "\n  Please modify the parameters N, minCellsPerBin, and/or alpha before attempting to resample"))
    }
    
  }
  
  # Determine and report back the average number of cells assigned to each unique bin
  # This is a useful metric for understanding how distinct the geometric sketch will be from random sampling
  # If the number of cells assigned to each unique bin is 1 then the geometric sketch is equivalent to random sampling
  cellsPerPlaid <- round(mean(table(Plaid)), 2)
  message(paste0(cellsPerPlaid, " cells per geometric bin."))
  
  sampledCells <- c() # Create the vector which sampled cells will be added to
  names(Plaid) <- row.names(X) # Ensure binIDs are associated with their respective cells
  
  # Sample only from the cells left in bins with more than the the minimum number of cells
  PlaidLeftover <- Plaid[Plaid %in% PlaidMinCells]
  
  sampledCells <- c() # Begin sampling cells
  while (length(sampledCells) < N) {
    
    # After each round of sampling ensure additional sampling is performed on unsampled cells only
    PlaidLeftover <- PlaidLeftover[!(names(PlaidLeftover) %in% sampledCells)]
    PlaidAddrsLeftover <- unique(PlaidLeftover)
    
    nCellsLeft <- N - (length(sampledCells)) # Calculate the number of cells left to sample
    if (nCellsLeft > length(PlaidAddrsLeftover)) { # If the number of cells left to sample is greater than the number of bins left, sample one cell from each bin
      sampleReamaining <- sapply(PlaidAddrsLeftover,
                                 FUN = function(PlaidAddr){names(sample(PlaidLeftover[PlaidLeftover == PlaidAddr], 1))})
    } else{ # Else sample the required number of cells one at a time from random remaining bins
      sampleReamaining <- sapply(sample(PlaidAddrsLeftover)[1:nCellsLeft],
                                 FUN = function(PlaidAddr){names(sample(PlaidLeftover[PlaidLeftover == PlaidAddr], 1))})
    }
    sampledCells <- c(sampledCells, sampleReamaining) # Add the newly sampled cells to the existing list of sampled cells
  }
  
  if (length(sampledCells) > N) { # In the event where the number of Bins was greater than N, randomly sample N cells from the already sampled cells
    sampledCells <- sample(sampledCells, N)
  }
  
  # Return sampled cellIDs and also cell binIDs if requested
  if (!returnBins) {
    return(sampledCells)
  } else {
    return(list(sampledCells = sampledCells, binIDs = Plaid))
  }
}
