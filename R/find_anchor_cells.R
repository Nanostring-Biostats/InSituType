#' Choose anchor cells
#' 
#' Finds cells with very good fits to the reference profiles, and saves these
#'  cells for use as "anchors" in the semi-supervised learning version of nbclust. 
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param bg Expected background
#' @param profiles Matrix of reference profiles holding mean expression of genes x cell types. 
#'  Input linear-scale expression, with genes in rows and cell types in columns.
#' @param size Negative binomial size parameter to be used in loglikelihood calculatoin
#' @param n_cells Up to this many cells will be taken as anchor points
#' @param min_cosine Cells must have at least this much cosine similarity to a fixed profile to be used as an anchor
#' @param min_scaled_llr Cells must have (log-likelihood ratio / totalcounts) above this threshold to be used as an anchor
#' @param insufficient_anchors_thresh Cell types that end up with fewer than this many anchors will be discarded. 
#' @return A vector holding anchor cell assignments (or NA) for each cell in the counts matrix
#' @importFrom lsa cosine
#' @export
find_anchor_cells <- function(counts, neg = NULL, bg = NULL, align_genes = TRUE,
                              profiles, size = 10, n_cells = 500, 
                              min_cosine = 0.3, min_scaled_llr = 0.01, 
                              insufficient_anchors_thresh = 20) {
  
  # infer bg if not provided: assume background is proportional to the scaling factor s
  if (is.null(bg) & is.null(neg)) {
    stop("Must provide either bg or neg")
  }
  if (is.null(bg)) {
    s <- Matrix::rowMeans(counts)
    bgmod <- stats::lm(neg ~ s - 1)
    bg <- bgmod$fitted
  }
  if (length(bg) == 1) {
    bg <- rep(bg, nrow(counts))
    names(bg) <- rownames(counts)
  }
  
  ### align genes in counts and fixed_profiles
  if (align_genes & !is.null(profiles)) {
    sharedgenes <- intersect(rownames(profiles), colnames(counts))
    lostgenes <- setdiff(colnames(counts), rownames(profiles))
    
    # subset:
    counts <- counts[, sharedgenes]
    profiles <- profiles[sharedgenes, ]
    
    # warn about genes being lost:
    if ((length(lostgenes) > 0) & length(lostgenes < 50)) {
      message(paste0("The following genes in the count data are missing from fixed_profiles and will be omitted from anchor selection: ",
                     paste0(lostgenes, collapse = ",")))
    }
    if (length(lostgenes) > 50) {
      message(paste0(length(lostgenes), " genes in the count data are missing from fixed_profiles and will be omitted from anchor selection"))
    }
  }
  
  # get cosine distances:
  cos <- matrix(NA, nrow(counts), ncol(profiles),
                dimnames = list(rownames(counts), colnames(profiles)))
  cos <- sapply(colnames(profiles), function(cell) {
    cos[, cell] <- apply(counts, 1, cosine, profiles[, cell])
  })
  
  # stats for which cells to get loglik on: 
  # get 3rd hightest cosines of each cell:
  cos3 <- apply(cos, 1, function(x){
    return(x[order(x, decreasing = TRUE)[3]])
  })
  # get cells with sufficient cosine:
  cells_with_high_cos <- apply(cos, 1, max) > min_cosine
  # get logliks (only when the cosine similarity is high enough to be worth considering):
  logliks <- sapply(colnames(profiles), function(cell) {
    templl <- cos[, cell] * NA
    usecells <- which((cos[, cell] >= pmin(0.75 * min_cosine, cos3)) & cells_with_high_cos)
    templl[usecells] <- Mstep(counts = counts[usecells, ], 
                                         means = profiles[, cell, drop = FALSE],
                                         freq = 1,
                                         bg = bg[usecells], 
                                         size = size, 
                                         digits = 3, return_loglik = T) 
    return(templl)
  })
  
  
  
  # choose anchors for each cell type:
  anchorslist <- lapply(colnames(profiles), function(cell) {
    
    usecells <- which(!is.na(logliks[, cell]))
    
    # get scaled log likelihood ratio:
    totcounts <- rowSums(counts[usecells, ])
    templlr <- logliks[usecells, cell] - apply(logliks[usecells, setdiff(colnames(logliks), cell)], 1, max, na.rm = TRUE)
    templlr <- templlr / totcounts
    tempcos <- cos[usecells, cell]
    # score based on llr and cosing
    score <- templlr * tempcos
    # require minimum values for each
    allowable <- which((templlr > min_scaled_llr) & (tempcos > min_cosine))
    useasanchors <- allowable[order(score[allowable], decreasing = TRUE)]
    useasanchors <- useasanchors[seq_len(min(n_cells, length(useasanchors)))]
    return(usecells[useasanchors])
  })
  names(anchorslist) <- colnames(profiles)
  anchors <- rep(NA, nrow(counts))
  names(anchors) <- rownames(counts)
  for (cell in names(anchorslist)) {
    anchors[anchorslist[[cell]]] <- cell
  }
  
  # anchor consolidation: identify and remove anchor cells that are poor fits to the mean anchor profile for the cell type:
  for (cell in setdiff(unique(anchors), NA)) {
    
    use <- (anchors == cell) & !is.na(anchors)
    # get centroid:
    if (!is.null(neg)) {
      mean_anchor_profile <- Estep(counts = counts[use, , drop = FALSE], clust = cell, neg = neg[use])
    } else {
      mean_anchor_profile <- Estep(counts = counts[use, , drop = FALSE], clust = cell, neg = bg[use])
    }
    
    # get anchors' cosine distances from centroid:
    newcos <- apply(counts[use, , drop = FALSE], 1, cosine, mean_anchor_profile)
    updated_anchors <- replace(anchors[use], (newcos < min_cosine), NA)
    anchors[names(updated_anchors)] <- updated_anchors
  }
  
  # remove all anchors from cells with too few total anchors:
  anchornums <- table(anchors)
  too_few_anchors <- names(anchornums)[anchornums <= insufficient_anchors_thresh]
  anchors[is.element(anchors, too_few_anchors)] <- NA
  if (length(too_few_anchors) > 0) {
    message(paste0("The following cell types had too few anchors and so are being removed from consideration: ",
                   paste0(too_few_anchors, collapse = ", ")))
  }

  return(anchors)
}
