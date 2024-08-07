#' Estimate the correct number of clusters using a subset of the data
#'
#' For a subset of the data, perform clustering under a range of cluster numbers.
#'  Report on loglikelihood vs. number of clusters, and suggest a best choice.
#' @param counts Counts matrix, cells * genes. 
#' @param neg Vector of mean negprobe counts per cell (default = "rna")
#' @param assay_type Assay type of RNA, protein 
#' @param bg Expected background
#' @param fixed_profiles Matrix of cluster profiles to hold unchanged throughout iterations.
#' @param fixed_sds Matrix of SDs expression of genes x cell types,to hold unchanged throughout iterations. Only for assay_type of protein
#' @param cohort Vector of cells' cohort assignments. 
#' @param init_clust Vector of initial cluster assignments.
#' @param n_clusts Vector giving a range of cluster numbers to consider.
#' @param max_iters Number of iterations in each clustering attempt. Recommended to choose
#'  a smaller number for a quicker, approximate clustering.
#' @param subset_size Number of cells to include in clustering.
#' @param align_genes Logical, for whether to align the genes in fixed_profiles with the colnames in count
#' @param plotresults Logical, for whether to plot the results.
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param pct_drop the decrease in percentage of cell types with a valid switchover to 
#'  another cell type compared to the last iteration. Default value: 1/10000. A valid 
#'  switchover is only applicable when a cell has changed the assigned cell type with its
#'  highest cell type probability increased by min_prob_increase. 
#' @param min_prob_increase the threshold of probability used to determine a valid cell 
#'  type switchover
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
#' @examples
#' data("mini_nsclc")
#' chooseClusterNumber(mini_nsclc$counts, Matrix::rowMeans(mini_nsclc$neg), assay_type="RNA",
#'  n_clust = 2:5)

chooseClusterNumber <-
  function(counts,
           neg,
           assay_type = c("rna", "protein"),
           bg = NULL,
           fixed_profiles = NULL,
           fixed_sds = NULL,
           cohort = NULL,
           init_clust = NULL,
           n_clusts = 2:12,
           max_iters = 10,
           subset_size = 1000,
           align_genes = TRUE,
           plotresults = FALSE,
           nb_size = 10,
           pct_drop = 0.005,
           min_prob_increase = 0.05,
           ...) {
    assay_type <- match.arg(tolower(assay_type), c("rna", "protein"))  

  # infer bg if not provided: assume background is proportional to the scaling factor s
  s <- rowSums(counts)
  if (is.null(bg)) {
    bgmod <- stats::lm(neg ~ s - 1)
    bg <- bgmod$fitted
  } 

  # subset the data:
  use <- sample(seq_len(nrow(counts)), subset_size)
  counts <- counts[use, ]
  s <- s[use]
  neg <- neg[use]
  bg <- bg[use]
  if (!is.null(init_clust)) {
    init_clust <- init_clust[use]
  }

  if (length(n_clusts) <= 0) {
    stop("n_clusts needs to be more than one value.")
  } else if (!all(sapply(n_clusts, function(x) x > 0 && x / as.integer(x) == 1))) {
    stop("n_clusts need to be a vector of positive integers.")
  }

  # align genes in fixed_profiles:
  if (align_genes && !is.null(fixed_profiles)) {
    sharedgenes <- intersect(rownames(fixed_profiles), colnames(counts))
    counts <- counts[, sharedgenes]
    fixed_profiles <- fixed_profiles[sharedgenes, ]
    fixed_sds <- fixed_sds[sharedgenes, ]
  }  
  # cluster under each value of n_clusts, and save loglik:
  totallogliks <- sapply(n_clusts, function(x) {
    
    # get init clust:
    tempinit <- rep(letters[seq_len(x)], each = ceiling(nrow(counts) / x))[
      seq_len(nrow(counts))]
   
    # run nbclust:
    message(sprintf("Clustering with n_clust = %s", x))
    tempclust <- nbclust(
      counts = counts, 
      neg = neg, 
      bg = bg, 
      fixed_profiles = fixed_profiles,
      fixed_sds = fixed_sds, 
      cohort = cohort,
      init_clust = tempinit,
      nb_size = nb_size,
      assay_type=assay_type,
      pct_drop = pct_drop,
      min_prob_increase = min_prob_increase,
      max_iters = max_iters)  

    # get the loglik of the clustering result:
    loglik_thisclust <- lldist(x = tempclust$profiles,
                               mat = counts,
                               xsd = tempclust$sds,
                               bg = bg,
                               size = nb_size,
                               assay_type = assay_type)

    total_loglik_this_clust <- sum(apply(loglik_thisclust, 1, max))
    return(total_loglik_this_clust)
  })

  # report goodness-of-fit
  n_parameters <- n_clusts * ncol(counts)
  aic <- n_parameters * 2 - 2 * totallogliks
  bic <- n_parameters * log(nrow(counts)) - 2 * totallogliks

  best_clust_number <- n_clusts[order(aic)[1]]

  if (plotresults) {
    original_par <- par()$mfrow
    graphics::par(mfrow = c(2, 1))
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
