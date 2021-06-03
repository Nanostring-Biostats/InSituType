if (FALSE) {
  counts = t(as.matrix(mini_tma@expression$rna$raw))
  s = Matrix::colSums(mini_tma@expression$rna$raw)
  bg = rep(0.1, length(s))
  celltype = mini_tma@cell_metadata$rna$slide
}


#' Mean background-subtracted expression per cluster
#' 
#' Estimates mean background-subtracted expression per cluster
#' @param counts Counts matrix, cells * genes.
#' @param s Vector of scaling factors for each cell, e.g. as defined by cell area. 
#' @param bg Expected background
#' @param celltype Vector of cell type assignments
#' @return A matrix of mean expression profiles
estimateCellTypeProfiles <- function(counts, s, bg, celltype) {
  
  # subtract background:
  scaledandsubtractedcounts <- sweep(sweep(counts, 1, bg, "-"), 1, s, "/")
  meanprofileslist <- by(scaledandsubtractedcounts, celltype, colMeans)
  mean_profiles <- pmax((matrix(unlist(meanprofileslist), ncol = length(meanprofileslist))), 0)
  colnames(mean_profiles) <- names(meanprofileslist)
  rownames(mean_profiles) <- colnames(scaledandsubtractedcounts)  
  return(mean_profiles)
}


#' Model platform effects between reference profiles and the current data
#' 
#' Estimates a scaling factor for each gene between the reference profiles and 
#'  the current data. Fits the model ______________________
#' @param counts Counts matrix, cells * genes.
#' @param s Vector of scaling factors for each cell, e.g. as defined by cell area. 
#' @param bg Expected background
#' @param celltype Vector of cell type assignments
#' @param fixed_profiles Matrix of expression profiles of pre-defined clusters,
#' @return A vector of each gene's calibration factor, the estimated ratio of efficiency in the data / efficiency in the fixed profile. 
estimateRefScalingFactors <- function(mean_profiles, fixed_profiles) {
  
  
  # Step 2: assemble the data frame for modeling:
  df <- data.frame(
    obs = as.vector(mean_profiles),
    ref = as.vector(fixed_profiles),
    gene = rep(rownames(mean_profiles), ncol(mean_profiles)),
    celltype = rep(colnames(mean_profiles), each = nrow(mean_profiles))
  )
  df$n <- table(celltype)[df$celltype]
  
  # calculate weights based on n and on counts (assuming poisson error with var = mean):
  df$wt <- sqrt(min(df$n, 1000) / pmax(df$obs, 1e-3))
  
  # run model:
  mod = lm(obs ~ ref*gene - 1, data = df, weights = df$wt) # this is not ideal 
  # better: run a separate model for each gene?
  # needed:
  # - adjust for cell type scaling (e.g. total RNA content)
  # - constrain so non-negative "ref" term
  
}






#' Estimate scaling factors between data and pre-specified reference profiles
#' 
#' a
#' @param counts Counts matrix, cells * genes.
#' @param s Vector of scaling factors for each cell, e.g. as defined by cell area. 
#' @param bg Expected background
#' @param celltype Vector of cell type assignments
#' @param fixed_profiles Matrix of expression profiles of pre-defined clusters,
#' @return A vector of each gene's calibration factor, the estimated ratio of efficiency in the data / efficiency in the fixed profile. 
estimateRefScalingFactors_stub <- function(counts, s, bg, celltype, fixed_profiles) {
  
  # For all genes, use the negbinom mle to get an optimal scaling factor.
  # the model: counts ~ NB(mean = r * s * x + b, theta = theta),
  #  where r = the gene's scaling factor, s = the vector of cells' scaling factors, x = the ref profile for the gene and b = background
  
  # problem: using this approach, we depend on having the cell profiles all on the right scale
  
  
  
  return(ref_scaling_factors)
}








#' Estimate scaling factors between data and pre-specified reference profiles
#' 
#' a
#' @param counts Counts matrix, cells * genes.
#' @param s Vector of scaling factors for each cell, e.g. as defined by cell area. 
#' @param bg Expected background
#' @param celltype Vector of cell type assignments
#' @param fixed_profiles Matrix of expression profiles of pre-defined clusters,
#' @return A vector of each gene's calibration factor, the estimated ratio of efficiency in the data / efficiency in the fixed profile. 
estimateRefScalingFactors_meanratio <- function(counts, s, bg, celltype, fixed_profiles) {
  
  # remove zeroes from fixed profiles:
  fixed_profiles = replace(fixed_profiles, fixed_profiles == 0, min(fixed_profiles[fixed_profiles > 0]))
  
  # init matrix of log fold-changes for each gene * cell type:
  meanfc <- matrix(NA, ncol(counts), ncol(fixed_profiles),
                   dimnames = list(colnames(counts), colnames(fixed_profiles)))
  
  # for all genes * cell types, get fold-change from reference to normalized data
  for (cell in colnames(meanfc)) {
    use = which(celltype == cell)
    if (length(use) > 50) {
      if (is.matrix(bg)) {
        meanexpr <- mean(sweep(pmax(counts[use, ] - bg[use, ], 0), 1, s[use], "/"))
      }
      if (length(bg) == nrow(counts)) {
        meanexpr <- colMeans(sweep(pmax(sweep(counts[use, ], 1, bg[use], "-"), 0), 1, s[use], "/"))
      }
      meanfc[, cell] <- meanexpr / fixed_profiles[, cell]
    }
    # should we scale the columns of meanfc to be on the same scale? 
    # -> this seems risky. Hard to prevent case where one column gets massively rescaled and skews the overall result. 
  }
  # rescale all of meanfc to avoid errors in small decimals:
  #meanfc = meanfc / quantile(meanfc, 0.9, na.rm = T)
  
  # get the weight for each cell type in each gene's average:
  # (note: the below calcultion is pretty arbitrary. It was chosen to advantage cell types with higher counts and with more cells)
  ncells <- table(celltype)[colnames(fixed_profiles)]
  wts <- sweep(sqrt(fixed_profiles), 2, sqrt(ncells), "*")
  wts[rowSums(wts) == 0, ] <- NA
  wts[meanfc == 0] <- 0
  wts <- sweep(wts, 1, rowSums(wts), "/")
  
  # for each gene, take a weighted average of log fold-changes, 
  ref_scaling_factors <- exp(rowSums(wts * (log(meanfc) - log(fixed_profiles)), na.rm = TRUE))
  return(ref_scaling_factors)
  
  # note 1: we get lots of meanFCs = 0, with log of -Inf. This will ruin calculations. 
  # x note 2: meanfc might hit scaling challenges - recenter everything to a reasonable place?
  # note 3: huge range of ref_scaling_factors. very long right tail. very implausible. either there's an error, or we need something more controlled like ebayes.
  # note 2/2: the huge range is probably driven in part by the ratios of very low-expressing genes. 
  # -> there should be some shrinkage applied. Perhaps this is a place to be bayesian.
  # -> Specify a prior based on CPA efficiencies, then update prior 
}
