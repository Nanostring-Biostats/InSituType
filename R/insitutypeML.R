#' Classify cells based on reference profiles
#' 
#' Supervised classification of cells. Each cell is assigned to the cell type 
#'  under which its observed expression profile is most likely. 
#' @param x Counts matrix (or dgCMatrix), cells * genes.
#'
#'   Alternatively, a \linkS4class{SingleCellExperiment} object containing such
#'   a matrix.
#' @param neg Vector of mean negprobe counts per cell. Can be provided 
#' @param bg Expected background
#' @param cohort Vector of cells' cohort memberships
#' @param reference_profiles Matrix of expression profiles of pre-defined clusters,
#'  e.g. from previous scRNA-seq. These profiles will not be updated by the EM algorithm.
#'  Colnames must all be included in the init_clust variable.
#' @param reference_sds Matrix of standard deviation profiles of pre-defined
#'   clusters. These SD profiles also will not be updated by the EM algorithm. 
#'   Columns must all be included in the init_clust variable. This parameter should
#'   be defined if assay_type is protein. Default is NULL. 
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param align_genes Logical, for whether to align the counts matrix and the reference_profiles by gene ID.
#' @param assay_type Assay type of RNA, protein (default = "rna") 
#' @param ... For the \linkS4class{SingleCellExperiment} method, additional
#'   arguments to pass to the ANY method.
#' @param assay.type A string specifying which assay values to use.
#' @return A list, with the following elements:
#' \enumerate{
#' \item clust: a vector given cells' cluster assignments
#' \item prob: a vector giving the confidence in each cell's cluster
#' \item profiles: Matrix of clusters' mean background-subtracted profiles
#' \item logliks: Matrix of cells' log-likelihoods under each cluster. Cells in rows, clusters in columns.
#' }
#'
#' @name insitutypeML
#' @examples
#' data("mini_nsclc")
#' data("ioprofiles")
#' sup <- insitutypeML(
#'  x = mini_nsclc$counts,
#'  neg = Matrix::rowMeans(mini_nsclc$neg),
#'  reference_profiles = ioprofiles,
#'  assay_type = "RNA")
#' table(sup$clust)
NULL

.insitutypeML <- function(x, neg = NULL, bg = NULL, cohort = NULL, 
                          reference_profiles, 
                          reference_sds=NULL, 
                          nb_size = 10, 
                          assay_type = c("rna", "protein"), 
                          align_genes = TRUE) {
  assay_type <- match.arg(tolower(assay_type), c("rna", "protein"))
  
  # get vector of expected background:
  bg <- estimateBackground(counts = x, neg = neg, bg = bg)
  
  # align genes:
  if (align_genes) {
    x <- alignGenes(counts = x, profiles = reference_profiles)
    reference_profiles <- reference_profiles[colnames(x), ]
    if (!is.null(reference_sds)) {
      reference_sds <- reference_sds[colnames(x), ]
    }
  }
  
  # prep cohort vector:
  if (is.null(cohort)) {
    cohort <- rep("all", length(bg))
  }
  
  logliks <- lldist(x = reference_profiles,
                    xsd = reference_sds,
                    mat = x,
                    bg =bg, 
                    size = nb_size,
                    assay_type=assay_type)
  

  # update logliks based on frequencies within cohorts:
  logliks <- update_logliks_with_cohort_freqs(logliks = logliks, 
                                              cohort = cohort, 
                                              minfreq = 1e-4, 
                                              nbaselinecells = 100) 
  if ("undefined" %in% colnames(logliks)) {
    logliks <- logliks[, -which(colnames(logliks) == "undefined")]
  }
  features <- intersect(rownames(reference_profiles), colnames(x))
  logliks <- cbind(logliks, ifelse(Matrix::rowSums(x[, features]) == 0, 0, -Inf))
  colnames(logliks)[ncol(logliks)] <- "undefined"
  
  # get remaining outputs
  clust <- colnames(logliks)[apply(logliks, 1, which.max)]
  names(clust) <- rownames(logliks)
  
  probs <- logliks2probs(logliks)
  prob <- apply(probs, 1, max)
  names(prob) <- names(clust)
  profiles_info <- Estep(counts=x, 
                         clust = clust,
                         neg = neg, 
                         assay_type=assay_type)
  
  profiles <- profiles_info$profiles
  sds <- profiles_info$sds
  
  # aligns profiles and logliks, removing lost clusters:
  logliks_from_lost_celltypes <- logliks[, !is.element(colnames(logliks), unique(clust)), drop = FALSE]
  logliks <- logliks[, is.element(colnames(logliks), clust), drop = FALSE]
  profiles <- profiles[, colnames(logliks), drop = FALSE]

  if(identical(tolower(assay_type), "rna")){
    sds <- NULL
  }
  
  out <- list(clust = clust,
             prob = prob,
             profiles = profiles,
             sds = sds,
             logliks = round(logliks, 4),
             logliks_from_lost_celltypes = round(logliks_from_lost_celltypes, 4))
  return(out)    
}

############################
# S4 method definitions 
############################

#' @export
#' @rdname insitutypeML
setGeneric("insitutypeML", function(x, ...) standardGeneric("insitutypeML"))

#' @export
#' @rdname insitutypeML
setMethod("insitutypeML", "ANY", .insitutypeML)

#' @export
#' @rdname insitutypeML
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment SingleCellExperiment
setMethod("insitutypeML", "SingleCellExperiment", function(x, ..., assay.type="counts") {
  .insitutypeML(t(assay(x, i=assay.type)), ...)
})
