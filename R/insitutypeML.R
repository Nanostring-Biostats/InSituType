#' Classify cells based on reference profiles
#' 
#' Supervised classification of cells. Each cell is assigned to the cell type 
#'  under which its observed expression profile is most likely. 
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell. Can be provided 
#' @param bg Expected background
#' @param fixed_profiles Matrix of expression profiles of pre-defined clusters,
#'  e.g. from previous scRNA-seq. These profiles will not be updated by the EM algorithm.
#'  Colnames must all be included in the init_clust variable.
#' @param align_genes Logical, for whether to align the counts matrix and the fixed_profiles by gene ID.
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param align_genes Logical, for whether to align the counts matrix and the fixed_profiles by gene ID.
#' @return A list, with the following elements:
#' \enumerate{
#' \item clust: a vector given cells' cluster assignments
#' \item probs: a matrix of probabilies of all cells (rows) belonging to all clusters (columns)
#' \item logliks: a matrix of each cell's log-likelihood under each cluster
#' }
#' @export
insitutypeML <- function(counts, neg = NULL, bg = NULL, fixed_profiles, nb_size = 10, align_genes = TRUE) {
  
  # get vector of expected background:
  if (is.null(bg) & is.null(neg)) {
    stop("Must provide either bg or neg")
  }
  # infer bg from neg if needed
  if (is.null(bg) & !is.null(neg)) {
      s <- Matrix::rowMeans(counts)
      bgmod <- stats::lm(neg ~ s - 1)
      bg <- bgmod$fitted
      names(bg) = rownames(counts)
  }
  # accept a single value of bg if input by user:
  if (length(bg) == 1) {
    bg <- rep(bg, nrow(counts))
    names(bg) <- rownames(counts)
  }
  
  # align genes:
  if (align_genes) {
    sharedgenes <- intersect(rownames(fixed_profiles), colnames(counts))
    lostgenes <- setdiff(colnames(counts), rownames(fixed_profiles))
    
    # subset:
    counts <- counts[, sharedgenes]
    fixed_profiles <- fixed_profiles[sharedgenes, ]
    
    # warn about genes being lost:
    if ((length(lostgenes) > 0) & length(lostgenes < 50)) {
      message(paste0("The following genes in the count data are missing from fixed_profiles and will be omitted from anchor selection: ",
                     paste0(lostgenes, collapse = ",")))
    }
    if (length(lostgenes) > 50) {
      message(paste0(length(lostgenes), " genes in the count data are missing from fixed_profiles and will be omitted from anchor selection"))
    }
  }
  
  # get logliks
  logliks <- apply(fixed_profiles, 2, function(ref) {
    lldist(x = ref,
           mat = counts,
           bg = bg,
           size = nb_size)
  })
  
  # get remaining outputs
  clust <- colnames(logliks)[apply(logliks, 1, which.max)]
  probs <- logliks2probs(logliks)
  out = list(clust = clust,
             probs = round(probs, 3),
             logliks = round(logliks, 3))
  return(out)    
  break()
}