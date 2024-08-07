#' @title Update cell typing results with spatial context or other alternative data
#' 
#' @description
#' Takes cell typing results, then updates it based on alternative data types, 
#' e.g. spatial context, morphology, or protein expression. Existing cell typing results are 
#' put into Insitutype's likelihood framework, which then can use alternative data
#' as a prior to be updated by the expression data to get a new posterior probability 
#' of cell type.
#' Performs this operation by 
#' \enumerate{
#' \item deriving cell type profiles using InSituType:::Estep(), 
#' \item assigning cells to "cohorts" (clusters) derived from their alternative data
#' \item  Inputing the output of steps (1) and (2) into InSituType::insitutype() to 
#'  re-calculate cell type. 
#' }
#' Paths for using alternative data in priority order (choose one; if multiple are input, only the most downstream option will be used):
#' \enumerate{
#' \item Input \code{xy} positions (and possibly \code{tissue}). Then cells will be clustered 
#'  into cohorts based on the expression pattern of their 50 nearest neighboring cells.
#' \item Input a matrix of alternative data (\code{altdata}) to be automatically clustered into cohorts. This supersedes 
#'  the altdata matrix derived from the \code{xy} argument.
#' \item Input your own \code{cohort} vector. This supersedes the above inputs. 
#' }
#' @param celltype Vector of cell type assignments to be updated
#' @param counts Counts matrix (or dgCMatrix), cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param cohort Vector of cells' cohort memberships. Output of a spatial clustering algorithm makes for good cohorts. 
#' @param altdata Matrix of cells' alternative data values
#' @param xy 2-column matrix of cells' xy positions. 
#' @param tissue Vector giving cells' tissue IDs. Used to separate tissue with overlapping xy coordinates.
#' @param nb_size The size parameter to assume for the NB distribution.
#' @param assay_type A string specifying which assay values to use.
#' @importFrom irlba irlba
#' @export
spatialUpdate <- function(celltype, counts, neg, 
                          cohort = NULL, altdata = NULL, xy = NULL, tissue = NULL,
                          nb_size = 10, assay_type = c("rna", "protein")) {
  
  assay_type <- match.arg(tolower(assay_type), c("rna", "protein"))
  
  ## check alternative data args:
  if(all(sapply(c(cohort, altdata, xy), is.null))) {
    stop("Must supply cohort, altdata or xy")
  }
  
  ## process alternative data, obtaining cohort vector:
  if (is.null(cohort)) {
    if (is.null(altdata)) {
      # make altdata from cells' neighborhoods:
      altdata <- getSpatialContext(counts = counts, xy = xy, tissue = tissue, 
                                   N = 50, rad = NULL, dim_reduce_to = 20) 
    }
    # cluster altdata to get cohort:
    cohort <- fastCohorting(mat = altdata, 
                            gaussian_transform = TRUE) 
  }
  
  ## derive reference profiles from initial cell type vector:
  profiles <- Estep(counts = counts, 
                    clust = celltype, 
                    neg = neg,
                    assay_type = assay_type)
  print(str(profiles))
  ## Run supervised cell typing with InSituType
  res <- insitutype(x = counts,
                    cohort = cohort,
                    neg = neg, 
                    reference_profiles = profiles$profiles,
                    reference_sds = profiles$sds,
                    n_clusts = 0,
                    update_reference_profiles = FALSE,
                    assay_type = assay_type)
  res$cohort <- cohort
  return(res)
}




