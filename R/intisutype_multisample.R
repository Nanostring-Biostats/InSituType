#' Unsupervised and semi-supervised clustering of single cell expression data.
#'
#' A wrapper for runinsitutype, to coordinate clustering across multiple tissues.
#' @param counts Counts matrix, cells * genes.
#' @param tissue Vector of cells' tissue IDs.
#' @param neg Vector of mean negprobe counts per cell
#' @param bg Expected background
#' @param anchors Vector giving "anchor" cell types, for use in semi-supervised clustering. 
#'  Vector elements will be mainly NA's (for non-anchored cells) and cell type names
#'  for cells to be held constant throughout iterations. 
#' @param n_clusts Number of clusters, in addition to any pre-specified cell types.
#'  Enter 0 to run purely supervised cell typing from fixed profiles. 
#'  Enter a range of integers to automatically select the optimal number of clusters. 
#' @param fixed_profiles Matrix of expression profiles of pre-defined clusters,
#'  e.g. from previous scRNA-seq. These profiles will not be updated by the EM algorithm.
#'  Colnames must all be included in the init_clust variable.
#' @param sketchingdata Optional matrix of data for use in non-random sampling via "sketching".
#'  If not provided, then the data's first 20 PCs will be used. 
#' @param align_genes Logical, for whether to align the counts matrix and the fixed_profiles by gene ID.
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param init_clust Vector of initial cluster assignments.
#' If NULL, initial assignments will be automatically inferred.
#' @param n_starts the number of iterations
#' @param n_benchmark_cells the number of cells for benchmarking
#' @param n_phase1 Subsample size for phase 1 (random starts)
#' @param n_phase2 Subsample size for phase 2 (refining in a larger subset)
#' @param n_phase3 Subsample size for phase 3 (getting final solution in a very large subset)
#' @param n_chooseclusternumber Subsample size for choosing an optimal number of clusters
#' @param pct_drop the decrease in percentage of cell types with a valid switchover to 
#'  another cell type compared to the last iteration. Default value: 1/10000. A valid 
#'  switchover is only applicable when a cell has changed the assigned cell type with its
#'  highest cell type probability increased by min_prob_increase. 
#' @param min_prob_increase the threshold of probability used to determine a valid cell 
#'  type switchover
#' @param max_iters Maximum number of iterations.
#' @param n_anchor_cells For semi-supervised learning. Maximum number of anchor cells to use for each cell type. 
#' @param min_anchor_cosine For semi-supervised learning. Cells must have at least this much cosine similarity to a fixed profile to be used as an anchor.
#' @param min_anchor_llr For semi-supervised learning. Cells must have (log-likelihood ratio / totalcounts) above this threshold to be used as an anchor
#' @param insufficient_anchors_thresh Cell types that end up with fewer than this many anchors after anchor selection will be discarded. 
#' @param anchor_replacement_thresh Threshold for calling whether a cluster has wandered away from its anchor cells. 
#'   Clusters whose anchors are reassigned at this rate or higher will be renamed, and only their anchor cells will be used to define a profile for the cluster.
#' @importFrom stats lm
#' @importFrom Matrix rowMeans
#' @importFrom Matrix colSums
#'
#' @return A list, with the following elements:
#' \enumerate{
#' \item clust: a vector given cells' cluster assignments
#' \item probs: a matrix of probabilies of all cells (rows) belonging to all clusters (columns)
#' \item profiles: a matrix of cluster-specific expression profiles
#' \item anchors: from semi-supervised clustering: a vector giving the identifies and cell types of anchor cells
#' }
#' @export


insitutype <- function(counts, tissue = NULL, neg, bg = NULL, 
                       anchors = NULL,
                       n_clusts,
                       fixed_profiles = NULL, 
                       sketchingdata = NULL,
                       align_genes = TRUE, nb_size = 10, 
                       init_clust = NULL, n_starts = 10, n_benchmark_cells = 50000,
                       n_phase1 = 5000, n_phase2 = 20000, n_phase3 = 100000,
                       n_chooseclusternumber = 2000,
                       pct_drop = 1/10000, min_prob_increase = 0.05, max_iters = 40,
                       n_anchor_cells = 500, min_anchor_cosine = 0.3, min_anchor_llr = 0.01, insufficient_anchors_thresh = 20,
                       anchor_replacement_thresh = 0.75) {
  
  # checks:
  if (is.null(tissue)) {
    tissue = rep("", nrow(counts))
  }
  if (length(tissue) != nrow(counts)) {
    stop("tissue vector should have length equal to nrow(counts).")
  }
  
  #### Using all tissues together, pick anchors for known cell types ------------------------
  ## get neg in condition 
  if (is.null(names(neg))) {
    names(neg) <- rownames(counts)
  }
  if (length(neg) != nrow(counts)) {
    stop("length of neg should equal nrows of counts.")
  }
  
  
  ### infer bg if not provided: assume background is proportional to the scaling factor s
  if (is.null(bg)) {
    s <- Matrix::rowMeans(counts)
    bgmod <- stats::lm(neg ~ s - 1)
    bg <- bgmod$fitted
    names(bg) = rownames(counts)
  }
  if (length(bg) == 1) {
    bg <- rep(bg, nrow(counts))
    names(bg) <- rownames(counts)
  }
  
  #### select anchor cells if not provided: ------------------------------
  if (is.null(anchors) & !is.null(fixed_profiles)) {
    message("automatically selecting anchor cells with the best fits to fixed profiles")
    anchors <- find_anchor_cells(counts = counts, 
                                 neg = NULL, 
                                 bg = bg, 
                                 profiles = fixed_profiles, 
                                 size = nb_size, 
                                 n_cells = n_anchor_cells, 
                                 min_cosine = min_anchor_cosine, 
                                 min_scaled_llr = min_anchor_llr,
                                 insufficient_anchors_thresh = insufficient_anchors_thresh) 
  }
  if (is.null(anchors) & any(n_clusts == 0)) {
    stop("No anchors were selected, and n_clusts = 0. The algorithm can't run under these conditions. 
         Solutions include: 1. make anchor selection more generous. 2. select anchors by hand. 3. increase n_clusts")
  }
  
  # test anchors are valid:
  anchorcellnames = NULL
  if (!is.null(anchors) & (length(anchors) != nrow(counts))) {
    stop("anchors must have length equal to the number of cells (row) in counts")
  }
  if (!is.null(anchors)) {
    names(anchors) <- rownames(counts)
    anchorcellnames <- names(anchors)[!is.na(anchors)]
  }
  
  
  
  #### For each tissue, run insitutype on its cells plus the anchors --------------------------
  
  resultslist <- list()
  for (tiss in unique(tissue)) {
    
    # cluster this tissue's cells along with all the anchors:
    use <- (tissue == tiss) | !is.na(anchors)
    
    tempres <- runinsitutype(
      counts = counts[use, ],
      neg = neg[use], 
      bg = bg[use], 
      anchors = anchors[use],
      n_clusts = n_clusts,
      fixed_profiles = fixed_profiles, 
      sketchingdata = sketchingdata,
      align_genes = align_genes, nb_size = nb_size, 
      init_clust = init_clust[use], 
      n_starts = n_starts, n_benchmark_cells = n_benchmark_cells,
      n_phase1 = n_phase1, n_phase2 = n_phase2, n_phase3 = n_phase3,
      n_chooseclusternumber = n_chooseclusternumber,
      pct_drop = pct_drop, min_prob_increase = min_prob_increase, max_iters = max_iters,
      n_anchor_cells = n_anchor_cells, min_anchor_cosine = min_anchor_cosine, 
      min_anchor_llr = min_anchor_llr, insufficient_anchors_thresh = insufficient_anchors_thresh,
      anchor_replacement_thresh = anchor_replacement_thresh
    )
    
    # rename new clusters by their tissue:
    
    unknownclusternames <- setdiff(tempres$clust, anchors)
    replacementnames <- paste0(tiss, "cl", seq_along(unknownclusternames))
    names(replacementnames) <- unknownclusternames
    
    for (oldname in unknownclusternames) {
      tempres$clust[tempres$clust == oldname] <- replacementnames[oldname]
      colnames(tempres$probs)[colnames(tempres$probs) == oldname] <- replacementnames[oldname]
      colnames(tempres$profiles)[colnames(tempres$profiles) == oldname] <- replacementnames[oldname]  
    }
    
    
    # keep only the results from this tissue's cells:
    keepcells <- rownames(counts)[tissue == tiss]
    tempres$clust <- tempres$clust[keepcells]
    tempres$probs <- tempres$probs[keepcells, ]
    tempres$anchors <- tempres$anchors[keepcells]
    
    resultslist[[tiss]] <- tempres
  }
  
  
  #### Combine results, labeling new clusters by tissue ------------------------------------
  
  # define empty output:
  uniqueclusts <- c()
  for (tiss in names(resultslist)) {
    uniqueclusts <- unique(c(uniqueclusts, resultslist[[tiss]]$clust))
  }
  uniqueclusts <- c(intersect(colnames(fixed_profiles), uniqueclust), 
                    setdiff(colnames(fixed_profiles), uniqueclust))
  
  out <- list()
  out$clust <- rep(NA, nrow(counts))
  names(out$clust) <- rownames(counts)
  out$probs <- matrix(NA, nrow(counts), length(uniqueclusts),
                      dimnames = list(rownames(counts), uniqueclusts))
  out$profiles <- matrix(NA, ncol(counts), length(uniqueclusts),
                         dimnames = list(colnames(counts), uniqueclusts))
  out$anchors <- anchors
  
  # fill in output:
  for (tiss in names(resultslist)) {
    tempcells <- names(resultslist[[tiss]]$clust)
    tempcelltypes <- colnames(resultslist[[tiss]]$probs)
    out$clust[tempcells] <- resultslist[[tiss]]$clust
    out$probs[tempcells, tempcelltypes, drop = FALSE] <- resultslist[[tiss]]$probs
  }
  
  # estimate profiles using all tissues:
  out$profiles <- Estep(counts, 
                        clust = out$clust,
                        neg = neg)
  
  return(out)
}