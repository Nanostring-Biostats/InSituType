
#' Update reference profiles
#'
#' Update reference profiles using pre-specified anchor cells, or if no anchors
#' are specified, by first choosing anchor cells
#' @param reference_profiles Matrix of reference profiles, genes * cell types
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param bg Expected background
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param anchors Vector giving "anchor" cell types, for use in semi-supervised
#'   clustering. Vector elements will be mainly NA's (for non-anchored cells)
#'   and cell type names for cells to be held constant throughout iterations.
#' @param n_anchor_cells For semi-supervised learning. Maximum number of anchor
#'   cells to use for each cell type.
#' @param min_anchor_cosine For semi-supervised learning. Cells must have at
#'   least this much cosine similarity to a fixed profile to be used as an
#'   anchor.
#' @param min_anchor_llr For semi-supervised learning. Cells must have
#'   (log-likelihood ratio / totalcounts) above this threshold to be used as an
#'   anchor
#' @return updated reference profiles
#' @export
#' @examples
#' data("mini_nsclc")
#' data("ioprofiles")
#' counts <- mini_nsclc$counts
#' astats <- get_anchor_stats(counts = mini_nsclc$counts,
#'  neg = Matrix::rowMeans(mini_nsclc$neg),
#'  profiles = ioprofiles)
#'
#' # estimate per-cell bg as a fraction of total counts:
#' negmean.per.totcount <- mean(rowMeans(mini_nsclc$neg)) / mean(rowSums(counts))
#' per.cell.bg <- rowSums(counts) * negmean.per.totcount
#'
#' # now choose anchors:
#' anchors <- choose_anchors_from_stats(counts = counts, 
#'                                     neg = mini_nsclc$negmean, 
#'                                     bg = per.cell.bg,
#'                                     anchorstats = astats, 
#'                                     # a very low value chosen for the mini
#'                                     # dataset. Typically hundreds of cells
#'                                     # would be better.
#'                                     n_cells = 50, 
#'                                     min_cosine = 0.4, 
#'                                     min_scaled_llr = 0.03, 
#'                                     insufficient_anchors_thresh = 5)
#'
#' # The next step is to use the anchors to update the reference profiles:
#'
#' updateReferenceProfiles(reference_profiles = ioprofiles, 
#'                        counts = mini_nsclc$counts, 
#'                        neg = mini_nsclc$neg, 
#'                        bg = per.cell.bg,
#'                        anchors = anchors) 
updateReferenceProfiles <-
  function(reference_profiles,
           counts,
           neg,
           bg = NULL,
           nb_size = 10,
           anchors = NULL,
           n_anchor_cells = 2000,
           min_anchor_cosine = 0.3,
           min_anchor_llr = 0.01) {
    
  
  
  ## step 1: if no anchors are provided, select them automatically:
  if (is.null(anchors)) {
    message("automatically selecting anchor cells with the best fits to fixed profiles")
    # align genes:
    sharedgenes <- intersect(colnames(counts), rownames(reference_profiles))
    anchors <- find_anchor_cells(counts = counts[, sharedgenes], 
                                 neg = neg, 
                                 bg = bg, 
                                 profiles = reference_profiles[sharedgenes, ], 
                                 size = nb_size, 
                                 n_cells = n_anchor_cells, 
                                 min_cosine = min_anchor_cosine, 
                                 min_scaled_llr = min_anchor_llr,
                                 insufficient_anchors_thresh = 0) 
  }
  
  if (is.null(anchors))  {
    stop("No anchors were selected. The algorithm can't run under these conditions. 
         Solutions include: 1. make anchor selection more generous. 2. select anchors by hand.")
  }
  
  # test anchors are valid:
  if (!is.null(anchors) && (length(anchors) != nrow(counts))) {
    stop("anchors must have length equal to the number of cells (row) in counts")
  }
  names(anchors) <- rownames(counts)
  
  ## step 2: use the anchors to update the reference profiles
  updated_profiles <- updateProfilesFromAnchors(counts = counts, 
                                                neg = neg, 
                                                anchors = anchors, 
                                                reference_profiles = reference_profiles, 
                                                align_genes = TRUE, nb_size = 10, max_rescaling = 5)
  out <- list(updated_profiles = updated_profiles,
              anchors = anchors)
  return(out)
}



#' Use anchor cells to update reference profiles, simply by taking the mean
#' profile of the anchors.
#'
#' Uses anchor cells to estimate platform effects / scaling factors to be
#' applied to the genes/rows of the reference profile matrix. Then uses Bayesian
#' math to update the individual elements on X.
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell. Can be provided
#' @param anchors Vector of anchor assignments
#' @param reference_profiles Matrix of expression profiles of pre-defined
#'   clusters, e.g. from previous scRNA-seq. These profiles will not be updated
#'   by the EM algorithm. Colnames must all be included in the init_clust
#'   variable.
#' @param align_genes Logical, for whether to align the counts matrix and the
#'   reference_profiles by gene ID.
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param max_rescaling Scaling factors will be truncated above by this value
#'   and below by its inverse (at 1/value and value)
#' @return \enumerate{ \item profiles: A profiles matrix with the rows rescaled
#' according to platform effects and individual elements updated further \item
#' scaling_factors: A vector of genes' scaling factors (what they were
#' multiplied by when updating the reference profiles). }
#' @export
updateProfilesFromAnchors <-
  function(counts,
           neg,
           anchors,
           reference_profiles,
           align_genes = TRUE,
           nb_size = 10,
           max_rescaling = 5) {
  use <- !is.na(anchors)
  updated_profiles <- Estep(counts = counts[use, ],
                            clust = anchors[use],
                            neg = neg[use])
  return(updated_profiles)
}
