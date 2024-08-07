#' Prepare bg data for other functions 
#' 
#' Process neg data or bg to get background for each cell
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param bg Expected background
#' @return A named vector for the estimated background of each cell

estimateBackground <- function(counts, neg, bg = NULL){
  # infer bg if not provided: assume background is proportional to the scaling factor s
  if (is.null(bg) && is.null(neg)) {
    stop("Must provide either bg or neg")
  }
  
  if (is.null(bg)) {
    ## get neg in condition 
    if (is.null(names(neg))) {
      names(neg) <- rownames(counts)
    }
    if (length(neg) != nrow(counts)) {
      stop("length of neg should equal nrows of counts.")
    }
    
    s <- Matrix::rowMeans(counts)
    bgmod <- stats::lm(neg ~ s - 1)
    bg <- bgmod$fitted
  }
  if (length(bg) == 1) {
    bg <- rep(bg, nrow(counts))
    names(bg) <- rownames(counts)
  }
  
  # overwrite if non-positive bg
  bg[bg <=0] <- min(1e-5, bg[bg>0])
  
  return(bg)
  
}


#' align genes in counts to profiles for other functions 
#' 
#' Process counts to have genes shared with profiles
#' @param counts Counts matrix, cells * genes.
#' @param profiles Matrix of reference profiles holding mean expression of genes x cell types. 
#'  Input linear-scale expression, with genes in rows and cell types in columns.
#' @return a cells * genes count matrix for shared genes only
alignGenes <- function(counts, profiles){
  sharedgenes <- intersect(rownames(profiles), colnames(counts))
  if (length(sharedgenes) < 10) {
    stop("Profiles have fewer than 10 genes in common with panel, use different profiles or re-run InSituType in unsupervised mode.")
  }
  lostgenes <- setdiff(colnames(counts), rownames(profiles))
  
  # subset:
  counts <- counts[, sharedgenes]
  
  # warn about genes being lost:
  if ((length(lostgenes) > 0) && length(lostgenes < 50)) {
    message(
      paste0(
        "The following genes in the count data are missing from fixed_profiles and will be omitted from downstream: ",
        paste0(lostgenes, collapse = ",")
      )
    )
  }
  if (length(lostgenes) > 50) {
    message(
      paste0(
        length(lostgenes),
        " genes in the count data are missing from fixed_profiles and will be omitted from downstream"
      )
    )
  }

  return(counts)
}


#' Get number of cores for parallelized operations
#'
#' @param percentCores percent of cores to use for parallelization [0-1]
#' @param minNotUsedCores minimum number of cores to leave for background processes
#' 
#' @return number of cores to use for mclapply
#' @export
numCores <- function(percentCores = 0.9, minNotUsedCores = 2) {
  if(percentCores > 1 & percentCores <= 0){
    stop("percentCores is not a valid number, must be between 0-1")
  }
  
  num_cores <- 1
  if (.Platform$OS.type == "unix") {
    if (is.null(getOption("mc.cores"))) {
      num_cores <- parallel::detectCores()
      if(num_cores <= minNotUsedCores){
        stop("minNotUsedCores must be fewer than available cores")
      }
      num_cores <- min(floor(num_cores*percentCores), num_cores-minNotUsedCores)
    } else {
      num_cores <- getOption("mc.cores") 
    }
  }
  return(num_cores)
}
