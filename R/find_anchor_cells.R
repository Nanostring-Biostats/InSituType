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
find_anchor_cells(counts, neg = NULL, bg = NULL, profiles, size = 10, n_cells = 500, min_cosine = 0.3, min_scaled_llr = 0.01) {
  
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
  
  # get logliks
  sharedgenes <- intersect(colnames(counts), rownames(profiles))
  logliks <- Mstep(counts = counts[, sharedgenes], 
                   means = profiles[sharedgenes, ],
                   bg = bg, 
                   size = size, 
                   digits = 2, return_loglik = T) 
  
  # get cosine distances:
  cos <- logliks * NA
  cos <- sapply(colnames(profiles), function(cell) {
    cos[, cell] <- apply(counts[, sharedgenes], 1, cosine, ioprofiles[sharedgenes, cell])
  })
  
  
  # choose anchors for each cell type:
  anchorslist <- lapply(colnames(profiles), function(cell) {
    # get scaled log likelihood ratio:
    totcounts <- rowSums(counts)
    llr <- logliks[, cell] - apply(logliks[, setdiff(colnames(logliks), cell)], 1, max)
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
  
  return(anchors)
}
