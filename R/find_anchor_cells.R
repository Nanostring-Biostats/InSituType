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
#' @return A vector holding anchor cell assignments (or NA) for each cell in the counts matrix
#' @importFrom lsa cosine
#' @export
find_anchor_cells <- function(counts, neg = NULL, bg = NULL, align_genes = TRUE, profiles, size = 10, n_cells = 500, min_cosine = 0.3, min_scaled_llr = 0.01) {
  
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
  
  # get logliks (only when the cosine similarity is high enough to be worth considering):
  cells_with_high_cos <- apply(cos, 1, max) > min_cosine
  
  logliks <- sapply(colnames(profiles), function(cell) {
    templl <- cos[, cell] * NA
    cells_to_consider <- which((cos[, cell] > 0.75 * min_cosine) & cells_with_high_cos)
    templl[cells_to_consider] <- Mstep(counts = counts[cells_to_consider, ], 
                                         means = profiles[, cell, drop = FALSE],
                                         freq = 1,
                                         bg = bg[cells_to_consider], 
                                         size = size, 
                                         digits = 3, return_loglik = T) 
    return(templl)
  })
  
  
  
  # choose anchors for each cell type:
  anchorslist <- lapply(colnames(profiles), function(cell) {
    # get scaled log likelihood ratio:
    totcounts <- rowSums(counts)
    llr <- logliks[, cell] - apply(logliks[, setdiff(colnames(logliks), cell)], 1, max, na.rm = TRUE)
    llr <- llr / totcounts
    # score based on llr and cosing
    score <- llr * cos[, cell]
    # require minimum values for each
    allowable <- which((llr > min_scaled_llr) & (cos[, cell] > min_cosine))
    useasanchors <- allowable[order(score[allowable], decreasing = TRUE)]
    useasanchors <- useasanchors[seq_len(min(n_cells, length(useasanchors)))]
    return(useasanchors)
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
      mean_anchor_profile <- Estep(counts = counts[use, ], clust = cell, neg = neg[use])
    } else {
      mean_anchor_profile <- Estep(counts = counts[use, ], clust = cell, neg = bg[use])
    }
    
    # get anchors' cosine distances from centroid:
    newcos <- apply(counts[use, ], 1, cosine, mean_anchor_profile)
    updated_anchors <- replace(anchors[use], (newcos < min_cosine), NA)
    anchors[names(updated_anchors)] <- updated_anchors
  }
  
  return(anchors)
}
