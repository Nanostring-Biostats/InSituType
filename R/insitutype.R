# NanoString Technologies, Inc.
# Software License Agreement for Non-Commercial Use
# By downloading, installing, accessing, modifying or otherwise making use of the Program (defined below), you agree to be bound by the terms and conditions of this Software License Agreement for Non-Commercial Use (this “License”).
# 1.	DEFINITIONS
# 1.1.	“Affiliate” means, with respect to an individual or entity, another individual or entity: (i) on whose behalf such individual or entity is acting, or (ii) that exercises control, is controlled by, or is under common control with such individual or entity. For the purposes of this definition, the term “control” means the right, whether by ownership, exercise of voting rights, contract, or otherwise, to direct the actions of an individual or entity.
# 1.2.	“Distribute” means to distribute, share, make available, or otherwise provide the Program or Modified Program, as applicable, or access thereto (including via a computer network) to any third party.
# 1.3.	“Licensor” means the individual or entity licensing the rights granted in this License.
# 1.4.	“Licensee” or “you” means the individual or entity receiving or exercising the rights granted under this License, provided that the individual or entity is not a NanoString Competitor.
# 1.5.	“Non-Commercial Use” means any use where profit or other commercial benefit is not a direct or indirect motive or intended result.
# 1.6.	“Modified Program” means a derivative work of, or a work that is based on, uses or incorporates, the Program (whether or not in combination with other works, materials or content).
# 1.7.	“NanoString” means NanoString Technologies, Inc.
# 1.8.	“NanoString Competitor” means any individual or entity that directly or indirectly competes with NanoString or any of NanoString’s Affiliates or whose Affiliate directly or indirectly competes with NanoString or any of NanoString’s Affiliates.
# 1.9.	“Program” means the copyrightable work of authorship, program, code, or software licensed under this License.
# 2.	LICENSE 
# 2.1.	Grant. Subject to the terms and conditions of this License, Licensor hereby grants to Licensee a worldwide, royalty-free, non-exclusive, revocable license to: (a) use, Distribute, and reproduce the Program, and (b) use, create, Distribute, and reproduce Modified Programs, in each case, solely for your internal, Non-Commercial Use. No rights are granted to NanoString Competitors.
# 2.2.	No Endorsement. Nothing in this License may be construed as permission to assert or imply that Licensor, NanoString, or other contributors to the Program sponsors, endorses, or is otherwise connected with the Licensee or the entity or institution that Licensee represents.
# 2.3.	Trademarks. Trademark rights are not licensed to you under this License.
# 2.4.	Grant of Patent License. Subject to the terms and conditions of this License, NanoString hereby grants to you a perpetual, worldwide, non-exclusive, no-charge, royalty-free, irrevocable (except as stated in this section) patent license to make, have made, use, import, and otherwise transfer the Program, where such license applies only to those patent claims licensable by NanoString that are necessarily infringed by Licensee alone or by combination of its modification(s) to the Program or Modified Program to which such modification(s) was submitted. If you institute patent litigation against any entity (including a cross-claim or counterclaim in a lawsuit) alleging that the Program, Modified Program, or a modification incorporated within the Program or a Modified Program constitutes direct or contributory patent infringement, then any patent licenses granted to you under this License for the Program or any such Modified Program shall terminate as of the date such litigation is filed.
# 3.	CONDITIONS TO THE RIGHT TO DISTRIBUTE
# 3.1.	Notices. If you Distribute the Program or a Modified Program in any form, you must also provide to the recipient:
# 3.1.1.	a copy of this License; and 
# 3.1.2.	for Modified Programs, prominent notices identifying the portions of the Modified Program that have been modified, stating that you have modified the Program.
# 3.2.	Attribution. Except as otherwise expressly permitted under this License, you must keep intact, and you may not modify or remove, any notices, disclaimers, or attributions included in or provided with the Program. In addition, you must also include a prominent hypertext link back to NanoString’s website at www.nanostring.com. 
# 3.3.	License. You may only Distribute the Program or the Modified Program under the terms of this License (or any later version, at your election). You may not offer or impose any additional or different terms or conditions that, or take any measures to, restrict the exercise of the rights granted under this License.
# 4.	NO REPRESENTATIONS OR WARRANTIES; LIMITATIONS OF LIABILITY
# 4.1.	Disclaimer. UNLESS OTHERWISE AGREED BY LICENSOR IN WRITING, TO THE FULLEST EXTENT PERMITTED BY APPLICABLE LAW, LICENSOR OFFERS THE PROGRAM AS-IS AND MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND WITH REGARD TO THE PROGRAM, WHETHER EXPRESS, IMPLIED, STATUTORY OR OTHERWISE, INCLUDING WITHOUT LIMITATION, WARRANTIES OF TITLE, MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NONINFRINGEMENT. THE LICENSOR DOES NOT REPRESENT OR WARRANT THAT THE PROGRAM WILL BE ERROR FREE AND DOES NOT PROMISE THAT ANY SUCH ERRORS WILL BE CORRECTED.
# SOME JURISDICTIONS DO NOT ALLOW FOR THE EXCLUSION OF IMPLIED WARRANTIES, SO THE FOREGOING MAY NOT APPLY TO YOU.
# 4.2.	Limitation of Liability. TO THE FULLEST EXTENT PERMITTED BY APPLICABLE LAW, IN NO EVENT WILL THE LICENSOR OR NANOSTRING BE LIABLE TO YOU UNDER ANY LEGAL THEORY FOR ANY DAMAGES OF ANY KIND, INCLUDING ANY SPECIAL, INCIDENTAL, CONSEQUENTIAL, PUNITIVE OR EXEMPLARY DAMAGES ARISING OUT OF OR RELATED TO THE PROGRAM OR USE THEREOF, EVEN IF LICENSOR OR NANOSTRING HAS BEEN ADVISED OF THE POSSIBILITY OR LIKELIHOOD OF SUCH DAMAGES.
# 5.	MISCELLANEOUS
# 5.1.	Right to Enforce. NanoString is an express third-party beneficiary of this License and will be entitled to enforce the provisions of this License as if it were a party hereto. 
# 5.2.	Waiver; Amendment. No term or provision hereof will be considered waived by the Licensor, and no breach excused by Licensor, unless such waiver or consent is in writing and signed by an authorized representative of Licensor.  The waiver by Licensor of, or consent by Licensor to, a breach of any provision of this License by the Licensee, will not constitute, operate or be construed as a waiver of, consent to, or excuse of any other or subsequent breach by Licensee.  This License may be amended or modified only by an agreement in writing signed by an authorized representative of each of Licensor and Licensee.

#' Run insitutype.
#'
#' A wrapper for nbclust, to manage subsampling and multiple random starts.
#' @param x Counts matrix (or dgCMatrix), cells * genes.
#'
#'   Alternatively, a \linkS4class{SingleCellExperiment} object containing such
#'   a matrix.
#' @param neg Vector of mean negprobe counts per cell
#' @param bg Expected background
#' @param assay_type Assay type of rna, protein (default = "rna")
#' @param anchors Vector giving "anchor" cell types, for use in semi-supervised
#'   clustering. Vector elements will be mainly NA's (for non-anchored cells)
#'   and cell type names for cells to be held constant throughout iterations.
#' @param cohort Vector of cells' cohort memberships
#' @param n_clusts Number of clusters, in addition to any pre-specified cell
#'   types. Enter 0 to run purely supervised cell typing from fixed profiles.
#'   Enter a range of integers to automatically select the optimal number of
#'   clusters.
#' @param reference_profiles Matrix of mean expression profiles of pre-defined
#'   clusters, e.g. from previous scRNA-seq. These profiles will not be updated
#'   by the EM algorithm. Columns must all be included in the init_clust
#'   variable.
#' @param reference_sds Matrix of standard deviation profiles of pre-defined
#'   clusters. These SD profiles also will not be updated by the EM algorithm. 
#'   Columns must all be included in the init_clust variable. This parameter should
#'   be defined if assay_type is protein. Default is NULL. 
#' @param update_reference_profiles Logical, for whether to use the data to
#'   update the reference profiles. Default and strong recommendation is TRUE.
#'   (However, if the reference profiles are from the same platform as the
#'   study, then FALSE could be better.)
#' @param sketchingdata Optional matrix of data for use in non-random sampling
#'   via "sketching". If not provided, then the data's first 20 PCs will be
#'   used.
#' @param align_genes Logical, for whether to align the counts matrix and the
#'   fixed_profiles by gene ID.
#' @param nb_size The size parameter to assume for the NB distribution. This 
#'    parameter is only for RNA.
#' @param init_clust Vector of initial cluster assignments. If NULL, initial
#'   assignments will be automatically inferred.
#' @param n_starts the number of iterations
#' @param n_benchmark_cells the number of cells for benchmarking
#' @param n_phase1 Subsample size for phase 1 (random starts)
#' @param n_phase2 Subsample size for phase 2 (refining in a larger subset)
#' @param n_phase3 Subsample size for phase 3 (getting final solution in a very
#'   large subset)
#' @param n_chooseclusternumber Subsample size for choosing an optimal number of
#'   clusters
#' @param pct_drop the decrease in percentage of cell types with a valid
#'   switchover to another cell type compared to the last iteration. Default
#'   value: 1/10000. A valid switchover is only applicable when a cell has
#'   changed the assigned cell type with its highest cell type probability
#'   increased by min_prob_increase.
#' @param min_prob_increase the threshold of probability used to determine a
#'   valid cell type switchover
#' @param max_iters Maximum number of iterations.
#' @param n_anchor_cells For semi-supervised learning. Maximum number of anchor
#'   cells to use for each cell type.
#' @param min_anchor_cosine For semi-supervised learning. Cells must have at
#'   least this much cosine similarity to a fixed profile to be used as an
#'   anchor.
#' @param min_anchor_llr For semi-supervised learning. Cells must have
#'   (log-likelihood ratio / totalcounts) above this threshold to be used as an
#'   anchor
#' @param insufficient_anchors_thresh Cell types that end up with fewer than
#'   this many anchors after anchor selection will be discarded.
#' @param refinement Logical, flag for further anchor refinement, used when update_reference_profiles = TRUE (default = FALSE)
#' @param rescale Logical, flag for platform effect correction, used when update_reference_profiles = TRUE (default = FALSE)
#' @param refit Logical, flag for fitting reference profiles to anchors, used when update_reference_profiles = TRUE (default = TRUE)
#' @param ... For the \linkS4class{SingleCellExperiment} method, additional
#'   arguments to pass to the ANY method.
#' @param assay.type A string specifying which assay values to use.
#' @importFrom stats lm
#' @importFrom Matrix rowMeans
#' @importFrom Matrix colSums
#' @export
#' @return A list, with the following elements: \enumerate{ \item clust: a
#'   vector given cells' cluster assignments \item prob: a vector giving the
#'   confidence in each cell's cluster \item logliks: Matrix of cells'
#'   log-likelihoods under each cluster. Cells in rows, clusters in columns.
#'   \item profiles: a matrix of cluster-specific expression profiles \item
#'   anchors: from semi-supervised clustering: a vector giving the identifies
#'   and cell types of anchor cells }
#'
#' @name insitutype
#' @examples 
#' data("mini_nsclc")
#' unsup <- insitutype(
#'  x = mini_nsclc$counts,
#'  neg = Matrix::rowMeans(mini_nsclc$neg),
#'  assay_type = "rna",
#'  n_clusts = 8,
#'  n_phase1 = 200,
#'  n_phase2 = 500,
#'  n_phase3 = 2000,
#'  n_starts = 1,
#'  max_iters = 5
#' ) # choosing inadvisably low numbers to speed the vignette; using the defaults in recommended.
#' table(unsup$clust)
NULL

#' @importFrom stats lm
#' @importFrom Matrix rowMeans
#' @importFrom Matrix colSums
#' @importFrom Matrix t
.insitutype <- function(x,
                        neg,
                        assay_type = c("rna", "protein"),
                        bg = NULL,
                        anchors = NULL,
                        cohort = NULL,
                        n_clusts,
                        reference_profiles = NULL,
                        reference_sds = NULL,
                        update_reference_profiles = TRUE,
                        sketchingdata = NULL,
                        align_genes = TRUE,
                        nb_size = 10,
                        init_clust = NULL,
                        n_starts = 10,
                        n_benchmark_cells = 10000,
                        n_phase1 = 10000,
                        n_phase2 = 20000,
                        n_phase3 = 100000,
                        n_chooseclusternumber = 2000,
                        pct_drop = 1 / 10000,
                        min_prob_increase = 0.05,
                        max_iters = 40,
                        n_anchor_cells = 2000,
                        min_anchor_cosine = 0.3,
                        min_anchor_llr = 0.03,
                        insufficient_anchors_thresh = 20,
                        refinement = FALSE, 
                        rescale = TRUE, 
                        refit = TRUE) {
  
  #### preliminaries ---------------------------

  assay_type <- match.arg(tolower(assay_type), c("rna", "protein"))
 
  temprows <- sample(seq_len(nrow(x)), min(c(1000, nrow(x))), replace = TRUE)
  tempcols <- sample(seq_len(ncol(x)), min(c(1000, nrow(x))), replace = TRUE)
  if ((assay_type == "rna") & any(abs(x[temprows, tempcols] - round(x[temprows, tempcols])) > 1e-4)) {
    warning("Non-integer elements of x input, and assay_type is set to rna. RNA mode Insitutype should use raw (i.e. integer) counts.")
  }
  
  if (any(rowSums(x) == 0)) {
    stop("Cells with 0 counts were found. Please remove.")
  }
  
  if (is.data.frame(reference_profiles)) {
    reference_profiles <- as.matrix(reference_profiles)
  }
  if (is.data.frame(reference_sds)) {
    reference_sds <- as.matrix(reference_sds)
  }
  
  # get vector of expected background:
  bg <- estimateBackground(counts = x, neg = neg, bg = bg)
  
  if (is.null(cohort)) {
    cohort <- rep("all", length(bg))
  }
  
  
  #### update reference profiles only for Supervised or Semi-Supervised clustering ----------------------------------
  fixed_profiles <- NULL
  fixed_sds <- NULL
  
  if (!is.null(reference_profiles)) { ## This is more for Supervised or Semi-Supervised case 
    if(is.null(reference_sds) && identical(tolower(assay_type), "protein")){
      stop("For protein data type, the reference SD profile should be provided!")
    }
    ## Update the profile matrix only for Semi-supervised or Supervised cases ##
    if (update_reference_profiles) {
      update_result <- updateReferenceProfiles(reference_profiles=reference_profiles,
                                               reference_sds=reference_sds,
                                               counts = x, 
                                               assay_type = assay_type,
                                               neg = neg,
                                               bg = bg,
                                               nb_size = nb_size,
                                               anchors = anchors,
                                               n_anchor_cells = n_anchor_cells, 
                                               min_anchor_cosine = min_anchor_cosine, 
                                               min_anchor_llr = min_anchor_llr,
                                               insufficient_anchors_thresh = insufficient_anchors_thresh, 
                                               refinement = refinement, 
                                               rescale = rescale, 
                                               refit = refit)
      fixed_profiles <- update_result$updated_profiles
      anchors <- update_result$anchors
      
      ## If assay_type==rna, this updated SDs are NULL
      fixed_sds <- update_result$updated_sds
      
    }else{
      fixed_profiles <- reference_profiles
      fixed_sds <- reference_sds
    }
  }
  
  # align the genes from fixed_profiles and counts
  if (align_genes && !is.null(fixed_profiles)) {
    x <- alignGenes(counts = x, profiles = fixed_profiles)
    fixed_profiles <- fixed_profiles[colnames(x), ]
    if(identical(tolower(assay_type), "protein")){
      fixed_sds <- fixed_sds[colnames(x), ]
    }
    if(identical(tolower(assay_type), "rna")){
      fixed_sds <- NULL
    }
  }

  
  #### set up subsetting: ---------------------------------
  # get data for subsetting if not already provided
  # (e.g., if PCA is the choice, then point to existing PCA results, and run PCA if not available
  if (!is.null(sketchingdata)) {
    # check that it's correct:
    if (nrow(sketchingdata) != nrow(x)) {
      warning("counts and sketchingdata have different numbers of row. Discarding sketchingdata.")
      sketchingdata <- NULL
    }
  }
  if (is.null(sketchingdata)) {
    sketchingdata <- prepDataForSketching(x, assay_type=assay_type)
  }
  n_phase1 <- min(n_phase1, nrow(x))
  n_phase2 <- min(n_phase2, nrow(x))
  n_phase3 <- min(n_phase3, nrow(x))
  n_benchmark_cells <- min(n_benchmark_cells, nrow(x))
  
  # define sketching "Plaids" (rough clusters) for subsampling:
  plaid <- geoSketch_get_plaid(X = sketchingdata, 
                               N = min(n_phase1, n_phase2, n_phase3, n_benchmark_cells),
                               alpha=0.1,
                               max_iter=200,
                               returnBins=FALSE,
                               minCellsPerBin = 1)
  
  #### choose cluster number: -----------------------------
  if (!is.null(init_clust)) {
    if (!is.null(n_clusts)) {
      message("init_clust was specified; this will overrule the n_clusts argument.")
    }
    n_clusts <- length(setdiff(unique(init_clust), colnames(fixed_profiles)))
  }
  if (is.null(n_clusts)) {
    n_clusts <- 10:15 + 5 * (is.null(fixed_profiles))
  }
  # get optimal number of clusters
  if (length(n_clusts) > 1) {
    
    message("Selecting optimal number of clusters from a range of ", min(n_clusts), " - ", max(n_clusts))

    chooseclusternumber_subset <- geoSketch_sample_from_plaids(Plaid = plaid, 
                                                               N = min(n_chooseclusternumber, nrow(x)))
    n_clusts <- chooseClusterNumber(
      counts = x[chooseclusternumber_subset, ], 
      neg = neg[chooseclusternumber_subset], 
      bg = bg[chooseclusternumber_subset], 
      fixed_profiles = reference_profiles,
      fixed_sds = reference_sds,
      init_clust = NULL, 
      n_clusts = n_clusts,
      max_iters = max(max_iters, 20),
      subset_size = length(chooseclusternumber_subset), 
      align_genes = TRUE, 
      plotresults = FALSE, 
      nb_size = nb_size, 
      assay_type=assay_type)$best_clust_number 
  }
  
  # This is for supervised case with or without reference profile update
  if(n_clusts==0){
    if(update_reference_profiles){
      profiles <- fixed_profiles
      sds <- fixed_sds
    }else{
      profiles <- reference_profiles
      sds <- reference_sds
    }
    
    if(is.null(reference_profiles)){
      
      stop("Reference profile should be provided for Supervised running.")
      
    }else{
      message(paste0("Supervised Case: Classifying all ", nrow(x), " cells with the user-defined reference profiles. "))
      
      out <- insitutypeML(x = x, 
                          neg = neg, 
                          bg = bg, 
                          reference_profiles = profiles, 
                          reference_sds = sds, 
                          cohort = cohort,
                          nb_size = nb_size,
                          assay_type=assay_type,
                          align_genes = TRUE) 
      
      ## Replace the cell types of anchor cells with the originally assigned anchor cells' cell types
      if(!is.null(anchors)){
        out$clust[!is.na(anchors)] <- anchors[names(out$clust[!is.na(anchors)])]
      }
      
      out$anchors <- anchors
      return(out)
      break
    }
  }
  
  #### phase 1: many random starts in small subsets -----------------------------
  if (!is.null(init_clust)) {
    message("init_clust was provided, so phase 1 - random starts in small subsets - will be skipped.")
    
    tempprofiles <- sapply(by(x[!is.na(init_clust), ], init_clust[!is.na(init_clust)], colMeans), cbind)
    rownames(tempprofiles) <- colnames(x)

    tempsds <- lapply(unique(init_clust[!is.na(init_clust)]), function(y){
      return(apply(x[!is.na(init_clust), ][init_clust[!is.na(init_clust)]==y,], 2, sd))
    })
    names(tempsds) <- unique(init_clust[!is.na(init_clust)])
    tempsds <- do.call("cbind", tempsds)
    
  } else {
    message(paste0("phase 1: random starts in ", n_phase1, " cell subsets"))
    
    # get a list in which each element is the cell IDs to be used in a subset
    random_start_subsets <- list()
    for (i in 1:n_starts) {
      random_start_subsets[[i]] <- geoSketch_sample_from_plaids(Plaid = plaid, 
                                                                N = min(n_phase1, nrow(x)))
    }
    
    # get a vector of cells IDs to be used in comparing the random starts:
    benchmarking_subset <- geoSketch_sample_from_plaids(Plaid = plaid, 
                                                        N = min(n_benchmark_cells, nrow(x)))
    
    # run nbclust from each of the random subsets, and save the profiles:
    profiles_from_random_starts <- list()
    sds_from_random_starts <- list()
    for (i in 1:n_starts) {
      
      cluster_name_pool <- c(letters, paste0(rep(letters, each = 26), rep(letters, 26)))
      tempinit <- rep(cluster_name_pool[seq_len(n_clusts)], each = ceiling(length(random_start_subsets[[i]]) / n_clusts))[
        seq_along(random_start_subsets[[i]])]

      tempNBclust <- nbclust(
        counts = x[random_start_subsets[[i]], ], 
        neg = neg[random_start_subsets[[i]]], 
        bg = bg[random_start_subsets[[i]]],
        fixed_profiles = fixed_profiles,
        fixed_sds = fixed_sds,
        cohort = cohort[random_start_subsets[[i]]],
        init_profiles = NULL,
        init_sds = NULL,
        init_clust = tempinit, 
        nb_size = nb_size,
        assay_type=assay_type,
        pct_drop = 1/500,
        min_prob_increase = min_prob_increase,
        max_iters = max(max_iters, 20),
      )
      profiles_from_random_starts[[i]] <- tempNBclust$profiles
      sds_from_random_starts[[i]] <- tempNBclust$sds
    }
    
    # find which profile matrix does best in the benchmarking subset:
    benchmarking_logliks <- c()
    for (i in 1:n_starts) {
      templogliks <- lldist(x = profiles_from_random_starts[[i]],
                            xsd = sds_from_random_starts[[i]],
                            mat = x[benchmarking_subset, ],
                            bg = bg[benchmarking_subset],
                            size = nb_size,
                            assay_type=assay_type)
      
      # take the sum of cells' best logliks:
      benchmarking_logliks[i] <- sum(apply(templogliks, 1, max))
    }
    best_start <- which.max(benchmarking_logliks)
    tempprofiles <- profiles_from_random_starts[[best_start]]
    
    if(identical(tolower(assay_type), "protein")){
      tempsds <- sds_from_random_starts[[best_start]]
    }
    
    if(identical(tolower(assay_type), "rna")){
      tempsds <- NULL
    }

    rm(profiles_from_random_starts)
    rm(sds_from_random_starts)
    rm(templogliks)
  }
  
  #### phase 2: -----------------------------------------------------------------
  message(paste0("phase 2: refining best random start in a ", n_phase2, " cell subset"))
  phase2_sample <- geoSketch_sample_from_plaids(Plaid = plaid, 
                                                N = min(n_phase2, nrow(x)))
  
  # get initial cell type assignments:
  temp_init_clust <- NULL
  if (!is.null(init_clust)) {
    temp_init_clust <- init_clust[phase2_sample]
    tempprofiles <- NULL
    tempsds <- NULL
  }
  
  # run nbclust, initialized with the cell type assignments derived from the previous phase's profiles
  clust2 <- nbclust(counts = x[phase2_sample, ], 
                    neg = neg[phase2_sample], 
                    bg = bg[phase2_sample],
                    fixed_profiles = fixed_profiles,
                    fixed_sds = fixed_sds,
                    cohort = cohort[phase2_sample],
                    init_profiles = tempprofiles, 
                    init_sds = tempsds, 
                    init_clust = temp_init_clust, 
                    nb_size = nb_size,
                    assay_type=assay_type,
                    pct_drop = 1/1000,
                    min_prob_increase = min_prob_increase,
                    max_iters = max_iters)
  tempprofiles <- clust2$profiles
  tempsds <- clust2$sds
  
  #### phase 3: -----------------------------------------------------------------
  message(paste0("phase 3: finalizing clusters in a ", n_phase3, " cell subset"))
  
  phase3_sample <- geoSketch_sample_from_plaids(Plaid = plaid, 
                                                N = min(n_phase3, nrow(x)))
  
  # run nbclust, initialized with the cell type assignments derived from the previous phase's profiles
  clust3 <- nbclust(counts = x[phase3_sample, ], 
                    neg = neg[phase3_sample], 
                    bg = bg[phase3_sample],
                    fixed_profiles = fixed_profiles,
                    fixed_sds = fixed_sds,
                    cohort = cohort[phase3_sample],
                    init_profiles = tempprofiles, 
                    init_sds = tempsds,
                    init_clust = NULL,  
                    nb_size = nb_size,
                    assay_type=assay_type,
                    pct_drop = pct_drop,
                    min_prob_increase = min_prob_increase,
                    max_iters = max_iters)
  profiles <- clust3$profiles
  sds <- clust3$sds
  
  
  #### phase 4: -----------------------------------------------------------------
  message(paste0("phase 4: classifying all ", nrow(x), " cells"))
  
  out <- insitutypeML(x = x, 
                      neg = neg, 
                      bg = bg, 
                      reference_profiles = profiles, 
                      reference_sds = sds,
                      cohort = cohort,
                      nb_size = nb_size, 
                      assay_type=assay_type,
                      align_genes = TRUE) 
  ## Replace the cell types of anchor cells with the originally assigned anchor cells' cell types
  if(!is.null(anchors)){
    out$clust[!is.na(anchors)] <- anchors[names(out$clust[!is.na(anchors)])]
  }
  out$anchors <- anchors
  
  return(out)
}

############################
# S4 method definitions 
############################

#' @export
#' @rdname insitutype
setGeneric("insitutype", function(x, ...) standardGeneric("insitutype"))

#' @export
#' @rdname insitutype
setMethod("insitutype", "ANY", .insitutype)

#' @export
#' @rdname insitutype
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment SingleCellExperiment
setMethod("insitutype", "SingleCellExperiment", function(x, ..., assay.type="counts") {
  .insitutype(t(assay(x, i=assay.type)), ...)
})