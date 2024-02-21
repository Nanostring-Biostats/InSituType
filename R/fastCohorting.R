#' Quickly split cells into cohorts 
#' 
#' Quickly split cells into cohorts using non-RNA data like spatial context and immunofluorescence values.
#' Rule of thumb: include any variables that might be informative for cell typing, 
#'  *except* variables you'll want to analyze later. For example, if you'll later
#'  perform differential expression as a function of spatial context, then it's 
#'  safer to exclude spatial context from the cell typing exercise (and therefore 
#'  from this function).
#' @param mat Matrix of variables to be used in cohorting, cells in rows, and variables in columns.
#'  Recommended to use < 20 variables. 
#' @param n_cohorts Number of clusters to divide cells into
#' @param gaussian_transform Whether to map each variable onto the quantiles of a normal distribution. 
#' @return A vector of cohort assignments. 
#' @export
#' @importFrom mclust Mclust
#' @importFrom mclust predict.Mclust
#' @importFrom mclust mclustBIC
#' @importFrom stats qnorm
#' @examples
#' data("mini_nsclc")
#' ## simulate immunofluorescence data: 
#' immunofluordata <- matrix(rpois(n = nrow(mini_nsclc$counts) * 4, lambda = 100), 
#'                           nrow(mini_nsclc$counts))
#' cohort <- fastCohorting(immunofluordata, gaussian_transform = TRUE)
#' table(cohort)
fastCohorting <- function(mat, n_cohorts = NULL, gaussian_transform = TRUE) {
  
  if (any(is.na(mat))) {
    stop("NA's detected in mat. fastCohorting needs complete data.")
  }

  # gaussian transform if called for:
  if (gaussian_transform) {
    for (i in seq_len(ncol(mat))) {
      mat[, i] <- qnorm(rank(mat[, i]) / (nrow(mat) + 1))
    }
  }
  
  # choose number of cohorts:
  if (is.null(n_cohorts)) {
    n_cohorts <- 3
    if (nrow(mat) > 10000) n_cohorts <- 10
    if (nrow(mat) > 50000) n_cohorts <- 25
    if (nrow(mat) > 100000) n_cohorts <- 50
    if (nrow(mat) > 200000) n_cohorts <- 100
  }
  
  # cluster in a subsample:
  sub <- sample(seq_len(nrow(mat)), min(20000, nrow(mat)))
  mc <- mclust::Mclust(data = mat[sub, ], G = n_cohorts, modelNames = "EEE")
  
  # classify all cells:
  cohort <- mclust::predict.Mclust(object = mc, newdata = mat)$classification
  names(cohort) <- rownames(mat)
  return(cohort)
}
