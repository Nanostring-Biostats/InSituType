#' Get anchor stats
#' 
#' Compute the statistics used in finding anchor cells.
#' Often the anchor cell selection process will involve some trial-and-error. 
#' This function performs the computationally-expensive steps that only need to 
#' happen once.
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param bg Expected background
#' @param align_genes Logical, for whether to align the columns of the counts matrix and the rows of
#'  the profiles matrix based on their names. 
#' @param profiles Matrix of reference profiles holding mean expression of genes x cell types. 
#'  Input linear-scale expression, with genes in rows and cell types in columns.
#' @param size Negative binomial size parameter to be used in likelihood calculation.
#' @param sds Matrix of reference profiles holding SDs expression of genes x cell types. 
#'  Input linear-scale expression, with genes in rows and cell types in columns. Only for assay_type of protein
#' @param assay_type Assay type of RNA, protein 
#' @param min_cosine Cells must have at least this much cosine similarity to a fixed profile to be used as an anchor.
#' @return A list with two elements: cos, the matrix of cosine distances;
#'  and llr, the matrix of log likelihood ratios of each cell under each cell type vs. the 2nd best cell type.
#' @importFrom lsa cosine
#' @export
#' @examples
#' data("ioprofiles")
#' data("mini_nsclc")
#' get_anchor_stats(counts = mini_nsclc$counts,
#'                  neg = Matrix::rowMeans(mini_nsclc$neg),
#'                  profiles = ioprofiles,
#'                  sds=NULL, 
#'                  assay_type = "RNA")
get_anchor_stats <- function(counts, neg = NULL, bg = NULL, align_genes = TRUE,
                             profiles, sds, size = 10, assay_type, 
                             min_cosine = 0.3) {
  
  # get vector of expected background:
  bg <- estimateBackground(counts = counts, neg = neg, bg = bg)
  
  ### align genes in counts and fixed_profiles
  if (align_genes && !is.null(profiles)) {
    counts <- alignGenes(counts = counts, profiles = profiles)
    profiles <- profiles[colnames(counts), ]
    sds <- sds[colnames(counts), ]
  }
  
  # get cosine distances:
  cos_numerator <- counts %*% as.matrix(profiles)
  counts2 <- as(counts, "dgCMatrix")
  counts2@x <- counts2@x^2
  rs <- sqrt(Matrix::rowSums(counts2))
  ps <- sqrt(Matrix::colSums(profiles^2))
  cos_denominator <- (t(t(rs))) %*% ps
  cos <- as.matrix(cos_numerator / cos_denominator)
  
  # stats for which cells to get loglik on: 
  # get 3rd hightest cosines of each cell:
  cos3 <- apply(cos, 1, function(x) {
    return(x[order(x, decreasing = TRUE)[3]])
  })
  
  # get cells with sufficient cosine:
  cells_with_high_cos <- apply(cos, 1, max) > min_cosine
  
  # get logliks (only when the cosine similarity is high enough to be worth considering):
  logliks <- sapply(colnames(profiles), function(cell) {
    templl <- cos[, cell] * NA
    usecells <- which((cos[, cell] >= pmin(0.75 * min_cosine, cos3)) & cells_with_high_cos)
    if (length(usecells) > 0) {
      templl[usecells] <- Mstep(counts = counts[usecells, ], 
                                means = profiles[, cell, drop = FALSE],
                                sds = sds[, cell, drop = FALSE],
                                cohort = rep("all", length(usecells)),
                                bg = bg[usecells], 
                                assay_type = assay_type,
                                size = size, 
                                digits = 3, return_loglik = TRUE) 
      # scale the logliks by total counts:
      templl[usecells] <- templl[usecells] / rowSums(counts[usecells,, drop = FALSE])
    }
    return(templl)
  })
  # convert loglik to LLR:
  getsecondbest <- function(x) {
    return(x[order(x, decreasing = TRUE)][2])
  }
  secondbest <- apply(logliks, 1, getsecondbest)
  llr <- suppressWarnings(sweep(logliks, 1, secondbest, "-"))
  
  out <- list(cos = cos, llr = llr)
  return(out)
}


#' Choose anchor cells given anchor stats
#'
#' Starting with cosine distances and log likelihood ratios, choose anchor
#' cells.
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param bg Expected background
#' @param anchorstats Output from get_anchor_stats. Must provide either this or
#'   both cos and llr matrices.
#' @param cos Matrix of cosine distances from reference profiles. Cells in rows,
#'   cell types in columns.
#' @param llr Matrix of log likelihood ratios from reference profiles. Cells in
#'   rows, cell types in columns.
#' @param n_cells Up to this many cells will be taken as anchor points
#' @param min_cosine Cells must have at least this much cosine similarity to a
#'   fixed profile to be used as an anchor
#' @param min_scaled_llr Cells must have (log-likelihood ratio / totalcounts)
#'   above this threshold to be used as an anchor
#' @param insufficient_anchors_thresh Cell types that end up with fewer than
#'   this many anchors will be discarded.
#' @param assay_type Assay type of RNA, protein 
#' 
#' @return A vector holding anchor cell assignments (or NA) for each cell in the
#'   counts matrix
#' @export
#' @examples 
#' data("ioprofiles")
#' data("mini_nsclc")
#' counts <- mini_nsclc$counts
#' astats <- get_anchor_stats(counts = counts,
#'                          neg = Matrix::rowMeans(mini_nsclc$neg),
#'                          sds=NULL, assay_type = "RNA",
#'                          profiles = ioprofiles)
#'
#' ## estimate per-cell bg as a fraction of total counts:
#' negmean.per.totcount <- mean(rowMeans(mini_nsclc$neg)) / mean(rowSums(counts))
#' per.cell.bg <- rowSums(counts) * negmean.per.totcount
#'
#' # now choose anchors:
#' choose_anchors_from_stats(counts = counts, 
#'                           neg = Matrix::rowMeans(mini_nsclc$neg),
#'                           bg = per.cell.bg,
#'                           anchorstats = astats, 
#'                           # a very low value chosen for the mini
#'                           # dataset. Typically hundreds of cells
#'                           # would be better.
#'                           n_cells = 50, 
#'                           min_cosine = 0.4, 
#'                           min_scaled_llr = 0.03, 
#'                           insufficient_anchors_thresh = 5,
#'                           assay_type="RNA")
choose_anchors_from_stats <-
  function(counts,
           neg = NULL,
           bg,
           anchorstats = NULL,
           cos = NULL,
           llr = NULL,
           n_cells = 500,
           min_cosine = 0.3,
           min_scaled_llr = 0.01,
           insufficient_anchors_thresh = 20,
           assay_type) {
    
    
    if (is.null(anchorstats) && (is.null(cos) || is.null(llr))) {
      stop("Must provide either anchorstats or both cos and llr matrices.")
    }

    # get input:
    if (!is.null(anchorstats)) {
      cos <- anchorstats$cos
      llr <- anchorstats$llr
    }
    
    # apply thresholds:
    cos <- cos * (cos > min_cosine)
    llr <- llr * (llr > min_scaled_llr)
    
    # choose anchors for each cell type:
    anchors <- rep(NA, nrow(cos))
    names(anchors) <- rownames(cos)
    for (cell in colnames(cos)) {
      score <- llr[, cell] * cos[, cell] * (llr[, cell] > min_scaled_llr) * (cos[, cell] > min_cosine)
      topn <- order(score, decreasing = TRUE)[seq_len(min(n_cells, sum(score > 0, na.rm = TRUE)))]
      rm(score)
      anchors[topn] <- cell
    }
    
    # anchor consolidation: identify and remove anchor cells that are poor fits to
    # the mean anchor profile for the cell type:
    for (cell in setdiff(unique(anchors), NA)) {
      
      use <- (anchors == cell) & !is.na(anchors)
      # get centroid:
      if (!is.null(neg)) {
        mean_anchor_profile <- Estep(counts = counts[use, , drop = FALSE], clust = cell, neg = neg[use], assay_type=assay_type)$profiles
      } else {
        mean_anchor_profile <- Estep(counts = counts[use, , drop = FALSE], clust = cell, neg = bg[use], assay_type=assay_type)$profiles
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
    
    if (length(setdiff(colnames(cos), unique(anchors))) > 0) {
      message(paste0("The following cell types had too few anchors and so are being removed from consideration: ",
                     paste0(setdiff(colnames(cos), unique(anchors)), collapse = ", ")))
    }
    
    if (all(is.na(anchors))) {
      warning("No anchor cells were selected - not enough cells met the selection criteria.")
      anchors <- NULL
    }
    return(anchors)  
  }





#' Choose anchor cells
#'
#' Finds cells with very good fits to the reference profiles, and saves these
#' cells for use as "anchors" in the semi-supervised learning version of
#' nbclust. The function would first pick anchor cell candidates through stats 
#' and then refine anchors based on umap projection. 
#' 
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param bg Expected background
#' @param align_genes Logical, for whether to align the columns of the counts
#'   matrix and the rows of the profiles matrix based on their names.
#' @param profiles Matrix of reference profiles holding mean expression of genes
#'   x cell types. Input linear-scale expression, with genes in rows and cell
#'   types in columns.
#' @param sds Matrix of reference profiles holding SDs expression of genes x cell types. 
#'  Input linear-scale expression, with genes in rows and cell types in columns. Only for assay_type of protein
#' @param size Negative binomial size parameter to be used in likelihood calculation. Only for assay_type of RNA
#' @param assay_type Assay type of RNA, protein 
#' @param n_cells Up to this many cells will be taken as anchor points
#' @param min_cosine Cells must have at least this much cosine similarity to a
#'   fixed profile to be used as an anchor
#' @param min_scaled_llr Cells must have (log-likelihood ratio / totalcounts)
#'   above this threshold to be used as an anchor
#' @param insufficient_anchors_thresh Cell types that end up with fewer than
#'   this many anchors will be discarded.
#' @param refinement flag to further refine the anchors via UMAP projection (default = FALSE)
#' @return A vector holding anchor cell assignments (or NA) for each cell in the
#'   counts matrix
#' @importFrom lsa cosine
#' @export
#' @examples 
#' data("ioprofiles")
#' data("mini_nsclc")
#' sharedgenes <- intersect(colnames(mini_nsclc$counts), rownames(ioprofiles))
#' find_anchor_cells(counts = mini_nsclc$counts[, sharedgenes], 
#'                   assay_type="RNA", 
#'                   sds=NULL,
#'                   neg = Matrix::rowMeans(mini_nsclc$neg),
#'                   profiles = ioprofiles)
find_anchor_cells <- function(counts, neg = NULL, bg = NULL, align_genes = TRUE,
                              profiles, sds, size = 10, assay_type,  n_cells = 500, 
                              min_cosine = 0.3, min_scaled_llr = 0.01, 
                              insufficient_anchors_thresh = 20,
                              refinement = FALSE) {
  
  if (align_genes && !is.null(profiles)) {
    counts <- alignGenes(counts = counts, profiles = profiles)
    profiles <- profiles[colnames(counts), ]
  }

  # get cos and llr stats:
  anchorstats <- get_anchor_stats(counts = counts,
                                  neg = neg,
                                  bg = bg, 
                                  align_genes = FALSE,
                                  profiles = profiles, 
                                  sds=sds, 
                                  size = size, 
                                  assay_type=assay_type, 
                                  min_cosine = min_cosine)  
  
  # select anchors based on stats:
  # double number for candidates if do further refinement
  anchors <- choose_anchors_from_stats(counts = counts, 
                                       neg = neg,
                                       bg = bg,
                                       anchorstats = anchorstats, 
                                       cos = NULL, 
                                       llr = NULL, 
                                       n_cells = ifelse(refinement, n_cells*2, n_cells), 
                                       min_cosine = min_cosine, 
                                       min_scaled_llr = min_scaled_llr, 
                                       insufficient_anchors_thresh = insufficient_anchors_thresh, 
                                       assay_type=assay_type) 
  
  if(refinement){
    # refine anchors via projection:
    anchors <- refineAnchors(counts = counts, 
                             neg = neg, bg = bg, 
                             align_genes =  FALSE,
                             profiles = profiles, 
                             anchor_candidates = anchors, 
                             nn_cells = n_cells,
                             insufficient_anchors_thresh = insufficient_anchors_thresh)
    
  }
  return(anchors)
}



#' Filter anchor candidates via projection of reference profiles to anchor-derived UMAP
#' 
#' Calculates expression UMAP model for anchor candidates, then projects reference 
#' profiles to the anchor-derived UMAP and select anchor candidates within top 
#' nearest neighbors of the projected reference profiles of same cell type in the 
#' UMAP as the final anchor cells. 
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param bg Expected background
#' @param align_genes Logical, for whether to align the columns of the counts matrix and the rows of
#'  the profiles matrix based on their names. 
#' @param profiles Matrix of reference profiles holding mean expression of genes x cell types. 
#'  Input linear-scale expression, with genes in rows and cell types in columns.
#' @param anchor_candidates Named vector of anchor candidates with cell_ID in name and corresponding cell type in values. 
#' @param nn_cells Number of top nearest neighbors to the projected reference profiles to be selected as final anchor cells. 
#' @param insufficient_anchors_thresh Cell types that end up with fewer than this many anchors will be discarded.
#' @return anchors, a named vector for the final anchor cells
#' @importFrom spatstat.geom ppp nncross
#' @importFrom uwot umap_transform
#' @export
refineAnchors <- function(counts, 
                          neg = NULL, 
                          bg = NULL, 
                          align_genes = TRUE,
                          profiles, 
                          anchor_candidates, 
                          nn_cells = 500, 
                          insufficient_anchors_thresh = 20) {
  # anchor candidates 
  cells_to_use <- names(anchor_candidates)[!is.na(anchor_candidates)]
  cells_to_use <- intersect(cells_to_use, rownames(counts))
  cts <- intersect(unique(anchor_candidates[cells_to_use]), colnames(profiles))
  
  if(length(cells_to_use)<20 || length(cts)<3){
    stop(sprintf("Only %d anchor candidates for %d cell types shared among `anchor_candidates`, `counts` and `profiles`. Must have at least 20 candidates for 3 cell types to use the projection based filtering.", 
                 length(cells_to_use), length(cts)))
  } else {
    message(sprintf("Start filtering on %d anchor candidates for %d cell types.", 
                    length(cells_to_use), length(cts)))
  }
  
  # infer bg if not provided: assume background is proportional to the scaling factor s
  bg <- estimateBackground(counts = counts, neg = neg, bg = bg)
  
  ### align genes in counts and fixed_profiles
  if (align_genes && !is.null(profiles)) {
    counts <- alignGenes(counts = counts, profiles = profiles)
    profiles <- profiles[colnames(counts), ]
  }
  
  # net expression profiles of anchor candidates, cell x gene 
  netExpr <- apply(sweep(as.matrix(counts[cells_to_use, , drop = F]), 1, bg[cells_to_use], "-"), 2, pmax, 0)
  netExpr <- netExpr[, Matrix::colSums(netExpr)>0, drop = F]
  
  sharedgenes <- intersect(rownames(profiles), colnames(netExpr))
  
  # get umap model for net expression of anchor candidates
  anc_umap <- uwot::umap(netExpr[, sharedgenes], metric = "cosine", ret_model = T)
  
  # projected ref in anc_umap
  projRef_umapcoord <- uwot::umap_transform(t(profiles)[, sharedgenes], anc_umap)
  
  # get top nearest neighbors for each projected reference within the anchors
  topNN <- as.matrix(spatstat.geom::nncross(
    # query point pattern for projected Ref
    X = spatstat.geom::ppp(x = projRef_umapcoord[, 1], 
                           y = projRef_umapcoord[, 2], 
                           range(projRef_umapcoord[, 1]), 
                           range(projRef_umapcoord[, 2]), 
                           marks = factor(rownames(projRef_umapcoord))), 
    # nearest neighbors in anchor data
    Y = spatstat.geom::ppp(x = anc_umap$embedding[, 1], 
                           y = anc_umap$embedding[, 2], 
                           range(anc_umap$embedding[, 1]), 
                           range(anc_umap$embedding[, 2]), 
                           marks = factor(anchor_candidates[rownames(anc_umap$embedding)])), 
    what = "which", k = 1:nn_cells))
  rownames(topNN) <- rownames(projRef_umapcoord)
  
  # get anchors within top Nearest Neighbors of consistent cell types
  anchors <- lapply(
    rownames(topNN),  
    function(ct){
      n500_cells <- rownames(anc_umap$embedding)[topNN[ct, ]]
      n500_cells <- intersect(n500_cells, 
                              names(anchor_candidates)[anchor_candidates == ct])
      ct_anchors <- rep(ct, length(n500_cells))
      names(ct_anchors) <- n500_cells
      return(ct_anchors)
    }
  )
  anchors <- do.call(c, anchors)
  
  message(sprintf("%d out of %d anchors are within top %d nearest neighbors of projected refProfiles of same cell types for %d cell types. ", 
                  length(anchors), length(cells_to_use), nn_cells, length(setdiff(unique(anchors), NA))))
  
  # remove all anchors from cells with too few total anchors:
  anchornums <- table(anchors)
  too_few_anchors <- names(anchornums)[anchornums <= insufficient_anchors_thresh]
  anchors[is.element(anchors, too_few_anchors)] <- NA
  
  if (length(setdiff(rownames(projRef_umapcoord), unique(anchors))) > 0) {
    message(paste0("The following cell types had too few anchors and so are being removed from consideration: ",
                   paste0(setdiff(rownames(projRef_umapcoord), unique(anchors)), collapse = ", ")))
  }
  
  if (all(is.na(anchors))) {
    warning("No anchor cells were selected - not enough cells met the selection criteria.")
    anchors <- NULL
  }  else {
    anchors <- setNames(anchors[names(anchor_candidates)], 
                        names(anchor_candidates))
  }
  
  return(anchors)
}


