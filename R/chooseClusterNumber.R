#' Estimate the correct number of clusters using a subset of the data
#'
#' For a subset of the data, perform clustering under a range of cluster numbers.
#'  Report on loglikelihood vs. number of clusters, and suggest a best choice.
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param fixed_profiles Matrix of expression profiles of pre-defined clusters,
#'  e.g. from previous scRNA-seq. These profiles will not be updated by the EM algorithm.
#'  Colnames must all be included in the init_clust variable.
#' @param init_clust Vector of initial cluster assignments.
#' @param n_clusts Vector giving a range of cluster numbers to consider.
#' @param n_iters Number of iterations in each clustering attempt. Recommended to choose
#'  a smaller number for a quicker, approximate clustering.
#' @param subset_size Number of cells to include in clustering.
#' @param align_genes Logical, for whether to align the genes in fixed_profiles with the colnames in count
#' @param plotresults Logical, for whether to plot the results.
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param ... Arguments passed to nbclust.
#' @export
#'
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom stats lm
#'
#' @return A list, with the following elements:
#' \itemize{
#'  \item
#' }
chooseClusterNumber <- function(counts, neg, bg = NULL, fixed_profiles = NULL,init_clust = NULL, n_clusts = 6:20,
                                n_iters = 10, subset_size = 1000, align_genes = TRUE, plotresults = FALSE, nb_size = 10, ...) {

  # infer bg if not provided: assume background is proportional to the scaling factor s
  if (is.null(bg)) {
    s <- rowSums(counts)
    bgmod <- stats::lm(neg ~ s - 1)
    bg <- bgmod$fitted
  }

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
  }

  # subset the data:
  set.seed(0)
  use <- sample(seq_len(nrow(counts)), subset_size)
  counts <- counts[use, ]
  s <- s[use]
  neg <- neg[use]
  bg <- bg[use]
  if (!is.null(init_clust)) {
    init_clust <- init_clust[use]
  }

  if (length(n_clusts) <=0 ){
    stop("n_clusts needs to be more than one value.")
  } else if ( !all( sapply(n_clusts, function(x) x>0 && x/as.integer(x) == 1) ) ){
    stop("n_clusts need to be a vector of positive integers.")
  }

  # cluster under each iteration, and save loglik:
  totallogliks <- sapply(n_clusts, function(x) {
    # run nbclust:
    message(sprintf("Fitting the model for n_clust = %s", x))
    tempclust <- nbclust(
      counts = counts, neg = neg, bg = bg,
      n_clusts = x, fixed_profiles = fixed_profiles,
      n_iters = n_iters)  # ,...

    # get the loglik of the clustering result:
    loglik_thisclust <- apply(tempclust$profiles, 2, function(ref) {
      lldist(x = ref,
             mat = counts,
             bg = bg,
             size = nb_size)
    })
    total_loglik_this_clust <- sum(apply(loglik_thisclust, 1, max))
  })

  # report goodness-of-fit
  n_parameters <- n_clusts * ncol(counts)
  aic <- n_parameters * 2 - 2 * totallogliks
  bic <- n_parameters * log(nrow(counts)) - 2 * totallogliks

  best_clust_number <- n_clusts[order(aic)[1]]

  if (plotresults) {
    original_par <- par()$mfrow
    graphics::par(mfrow = c(2,1))
    graphics::plot(n_clusts, totallogliks, xlab = "Number of clusters", ylab = "Log-likelihood")
    graphics::lines(n_clusts, totallogliks)
    graphics::plot(n_clusts, aic, xlab = "Number of clusters", ylab = "AIC")
    graphics::lines(n_clusts, aic)
    par(mfrow = original_par)
  }

  out <- list(best_clust_number = best_clust_number,
              n_clusts = n_clusts,
              loglik = totallogliks,
              aic = aic,
              bic = bic)
  return(out)


}

