

#' Estimate platform effects / gene scaling factors
#'
#' Estimates a scaling factor for each gene to bring the fixed_profiles onto the scale of the data.
#'  Models the log ratio between platforms as a function of gene, cell line, and
#'  reference expression.
#' @param mean_profiles Matrix of mean background-subtracted expression of genes * clusters.
#' @param fixed_profiles Matrix of pre-specified reference profiles.
#'
#' @importFrom stats quantile
#' @importFrom stats lm
#'
#' @return A vector of scaling factors used to multiply the rows of fixed_profiles
estimate_platform_effects <- function(mean_profiles, fixed_profiles) {

  if(!identical(dim(mean_profiles), dim(fixed_profiles))){
    stop("Please ensure the dimensions of mean_profiles and fixed_profiles are the same.")
  }

  # make data frame for model:
  df <- data.frame(obs = as.vector(mean_profiles),
                  ref = as.vector(fixed_profiles),
                  gene = rep(rownames(mean_profiles), ncol(mean_profiles)),
                  cell = rep(colnames(mean_profiles), each = nrow(mean_profiles)))


  ## find columns of mean_profiles with insufficient info, and remove:
  # flag all-0 columns:
  all_zero_columns <- colnames(mean_profiles)[colSums(mean_profiles) == 0]
  # flag columns with too few non-zero values:
  prop_nonzero <- colMeans(mean_profiles > 0)
  poor_data_columns <- colnames(mean_profiles)[prop_nonzero < 0.1 * mean(prop_nonzero)]
  # remove:
  df <- df[!is.element(df$cell, c(all_zero_columns, poor_data_columns)), ]

  # remove rows where ref == 0:
  df <- df[df$ref > 0, ]

  # lower-threshold the observed means so 0 values don't go to -Inf:
  
  df$obs <- pmax(df$obs, stats::quantile(df$obs[df$obs>0], 0.01, na.rm = TRUE))

  # calc logratio
  df$logratio <- log(df$obs / df$ref)

  # fit the model:
  mod <- stats::lm(logratio ~ 0 + gene + cell, data = df, weights = log(ref) - min(log(ref)))
  coefs <- mod$coefficients[grepl("gene", names(mod$coefficients))]

  # add coefs for genes that got dropped due to all 0s in df:
  names(coefs) <- gsub("gene", "", names(coefs))
  lostgenes <- setdiff(rownames(mean_profiles), names(coefs))
  coefs[lostgenes] <- mean(coefs)
  coefs <- coefs[rownames(mean_profiles)]
  return(exp(coefs))
}




#' Mean background-subtracted expression per cluster
#'
#' Estimates mean background-subtracted expression per cluster
#' @param counts Counts matrix, cells * genes.
#' @param s Vector of scaling factors for each cell, e.g. as defined by cell area.
#' @param bg Expected background
#' @param celltype Vector of cell type assignments
#' @return A matrix of mean expression profiles
#estimateCellTypeProfiles <- function(counts, s, bg, celltype) {
#
#  # subtract background:
#  scaledandsubtractedcounts <- sweep(sweep(counts, 1, bg, "-"), 1, s, "/")
#  meanprofileslist <- by(scaledandsubtractedcounts, celltype, colMeans)
#  mean_profiles <- pmax((matrix(unlist(meanprofileslist), ncol = length(meanprofileslist))), 0)
#  colnames(mean_profiles) <- names(meanprofileslist)
#  rownames(mean_profiles) <- colnames(scaledandsubtractedcounts)
#  return(mean_profiles)
#}




