
#' Update reference profiles
#' 
#' Update reference profiles using pre-specified anchor cells, or if no anchors are specified, by first choosing anchor cells
updateReferenceProfiles <- function(reference_profiles, counts, neg, 
                                      anchors = NULL, n_anchor_cells = 2000, min_anchor_cosine = 0.3, min_anchor_llr = 0.01) {
  
  
  ## step 1: if no anchors are provided, select them automatically:
  if (is.null(anchors)) {
    message("automatically selecting anchor cells with the best fits to fixed profiles")
    # align genes:
    sharedgenes <- intersect(colnames(counts), rownames(reference_profiles))
    anchors <- find_anchor_cells(counts = counts[, sharedgenes], 
                                 neg = NULL, 
                                 bg = bg, 
                                 profiles = reference_profiles[sharedgenes, ], 
                                 size = nb_size, 
                                 n_cells = n_anchor_cells, 
                                 min_cosine = min_anchor_cosine, 
                                 min_scaled_llr = min_anchor_llr,
                                 insufficient_anchors_thresh = insufficient_anchors_thresh) 
  }
  
  if (is.null(anchors))  {
    stop("No anchors were selected. The algorithm can't run under these conditions. 
         Solutions include: 1. make anchor selection more generous. 2. select anchors by hand.")
  }
  
  # test anchors are valid:
  if (!is.null(anchors) & (length(anchors) != nrow(counts))) {
    stop("anchors must have length equal to the number of cells (row) in counts")
  }
  names(anchors) <- rownames(counts)
  
  ## step 2: use the anchors to update the reference profiles
  updated_profiles <- updateProfilesFromAnchors(counts = counts, 
                                                neg = neg, 
                                                anchors = anchors, 
                                                reference_profiles = reference_profiles, 
                                                align_genes = TRUE, nb_size = 10, max_rescaling = 5)
  return(updated_profiles)
}

#' Rescale rows of reference profiles
#' 
#' Estimates platform effects and rescales the reference matrix accordingly
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell. Can be provided 
#' @param reference_profiles Matrix of expression profiles of pre-defined clusters,
#'  e.g. from previous scRNA-seq. These profiles will not be updated by the EM algorithm.
#'  Colnames must all be included in the init_clust variable.
#' @param align_genes Logical, for whether to align the counts matrix and the reference_profiles by gene ID.
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param max_rescaling Scaling factors will be truncated above by this value and below by its inverse (at 1/value and value)
#' @return 
#' \enumerate{
#' \item profiles: A profiles matrix with the rows rescaled according to platform effects
#' \item scaling_factors: A vector of genes' scaling factors (what they were multiplied by when updating the reference profiles). 
#' }
#' @export
#' @importFrom SpatialDecon spatialdecon
rescaleProfiles <- function(counts, neg, reference_profiles, align_genes = TRUE, nb_size = 10, max_rescaling = 5) {
  
  # align genes:
  if (align_genes) {
    sharedgenes <- intersect(rownames(reference_profiles), colnames(counts))
    lostgenes <- setdiff(colnames(counts), rownames(reference_profiles))
    
    # subset:
    counts <- counts[, sharedgenes]
    reference_profiles <- reference_profiles[sharedgenes, ]
    
    # warn about genes being lost:
    if ((length(lostgenes) > 0) & length(lostgenes < 50)) {
      message(paste0("The following genes in the count data are missing from reference_profiles and will be omitted from anchor selection: ",
                     paste0(lostgenes, collapse = ",")))
    }
    if (length(lostgenes) > 50) {
      message(paste0(length(lostgenes), " genes in the count data are missing from reference_profiles and will be omitted from anchor selection"))
    }
  }
  
  # get background-subtracted bulk profile:
  totcounts <- matrix(colSums(counts), ncol(counts))
  rownames(totcounts) <- colnames(counts)
  totcountsbgsub <- pmax(totcounts - sum(neg), 0)
  
  # deconvolve the bulk profile with spatialdecon:
  res <- SpatialDecon::spatialdecon(norm = cbind(totcounts, totcounts), 
                      bg = sum(neg), 
                      X = reference_profiles, 
                      resid_thresh = Inf)
  log2resids <- res$resids[, 1]
  log2resids <- log2resids - mean(log2resids)
  scaling_factors <- 2^log2resids
  scaling_factors <- pmax(pmin(scaling_factors, max_rescaling), max_rescaling^-1)
  rescaledprofiles <- sweep(reference_profiles, 1, scaling_factors, "*")
  
  out = list(profiles = rescaledprofiles,
             scaling_factors = scaling_factors)
  return(out)
}
  



#' Use anchor cells to update reference profiles
#' 
#' Uses anchor cells to estimate platform effects / scaling factors to be applied to
#'  the genes/rows of the reference profile matrix. Then uses Bayesian math to
#'  update the individual elements on X. 
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell. Can be provided 
#' @param anchors Vector of anchor assignments
#' @param reference_profiles Matrix of expression profiles of pre-defined clusters,
#'  e.g. from previous scRNA-seq. These profiles will not be updated by the EM algorithm.
#'  Colnames must all be included in the init_clust variable.
#' @param align_genes Logical, for whether to align the counts matrix and the reference_profiles by gene ID.
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param max_rescaling Scaling factors will be truncated above by this value and below by its inverse (at 1/value and value)
#' @return 
#' \enumerate{
#' \item profiles: A profiles matrix with the rows rescaled according to platform effects and individual elements updated further
#' \item scaling_factors: A vector of genes' scaling factors (what they were multiplied by when updating the reference profiles). 
#' }
#' @export
updateProfilesFromAnchors <- function(counts, neg, anchors, reference_profiles, align_genes = TRUE, nb_size = 10, max_rescaling = 5) {
  
  # get total expression and negs per anchor cell type:
  sums <- sapply(by(counts, anchors, colSums), cbind)
  rownames(sums) <- colnames(counts)
  temp <- by(neg, anchors, sum)
  negsums <- as.vector(temp)
  names(negsums) <- names(temp)
  negsums <- negsums[colnames(sums)]
  
  ### build data frame for estimating row/gene scaling factors:
  celltypes <- intersect(colnames(reference_profiles), unique(anchors))
  genes <- intersect(rownames(sums), rownames(reference_profiles))
  df <- data.frame(
    gene = rep(genes, length(celltypes)),
    celltype = rep(celltypes, each = lenght(genes)),
    ref = as.vector(reference_profiles[genes, celltypes]),
    rnasum = as.vector(sums[genes, celltypes]),
    negsum = rep(negsums, each = length(genes))
  )
  df$bgsub <- pmax(df$rnasum - df$negsum, 0)
  # throw out rows with values too low to be useful in estimating a log-ratio:
  
  ### calculate weights based on precision of anchors' total counts:
  # get sd of log(totrna - totneg) - log(ref), based on the delta method
  df$sd <- sqrt(sum(c((df$rnasum-df$negsum)^-2, (df$rnasum-df$negsum)^-2, df$ref^-2) * c(df$rnasum, df$negsum, df$ref)))
  
  #get_sd_of_log_bgsub_totcounts <- function(i) {
  #  y <- df$rnasum[i]
  #  n <- df$negsum[i]
  #  r <- df$ref[i]
  #  out = sqrt(t(c((y-n)^-1, -(y-n)^-1, -r^-1)) %*% diag(c(y, n, r)) %*% c((y-n)^-1, -(y-n)^-1, -r^-1))
  #  
  #  sqrt(sum(c((y-n)^-2, (y-n)^-2, r^-2) * c(y,n,r)))
  #  #out = sqrt(t(c((y-n-r)^-1, -(y-n-r)^-1, -(y-n-r)^-1)) %*% diag(c(sqrt(y), sqrt(n), sqrt(r))) %*% (c((y-n-r)^-1, -(y-n-r)^-1, -(y-n-r)^-1)))
  #  return()
  #}
  #df$sd = sapply(seq_len(nrow(df)), get_sd_of_log_bgsub_totcounts)
  df$weight = 1 / df$sd
  
  ### model platform effects:
  # remove rows with values too low to be useful:
  bgsub_toolow <- df$bgsub < 50
  ref_toolow <- df$ref < quantile(df$ref[df$ref > 0], 0.2)
  df <- df[!bgsub_toolow & !ref_toolow, ]
  mod <- lm(log(bgsub) ~ gene + celltype + offset(log(ref)) - 1, weights = df$weight, data = df)
  coefs <- mod$coefficients[substr(names(mod$coefficients), 1, 4) == "gene"]
  stderrs <- summary(mod)$coefficients[names(coefs), 2]
  # get a shrinkage estimate of coefs, assuming that coefs have a distribution with sd = log(1.25):
  priorvar <- log(1.25)^2
  coefs <- (coefs - mean(coefs)) * priorvar / (stderrs^2 + priorvar)
  # fill in missing genes:
  names(coefs) <- substr(names(coefs), 5, nchar(names(coefs)))
  coefs[is.na(coefs)] <- mean(coefs, na.rm = TRUE)
  coefs[setdiff(genes, names(coefs))] <- mean(coefs)
  
  # convert to log2-scale:
  log2_scaling <- log2(exp(1)) * coefs
  
  
  
}