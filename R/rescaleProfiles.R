
#' Update reference profiles
#'
#' Update reference profiles using pre-specified anchor cells, or if no anchors
#' are specified, by first choosing anchor cells. Option to return reference 
#' profiles rescaled for platform effect and/or to return further refitted profiles 
#' based on the observed profiles of anchor cells.
#' 
#' @param reference_profiles Matrix of reference mean profiles, genes * cell types
#' are specified, by first choosing anchor cells.
#' @param reference_sds Matrix of standard deviation profiles, genes * cell types. Only for assay_type of protein.
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param assay_type Assay type of RNA, protein (default = "rna")
#' @param bg Expected background
#' @param nb_size The size parameter to assume for the NB distribution. Only for assay_type of RNA
#' @param anchors named vector giving "anchor" cell types with cell_id in names, 
#' for use in semi-supervised clustering. Vector elements will be mainly NA's 
#' (for non-anchored cells) and cell type names for cells to be held constant 
#' throughout iterations.
#' @param n_anchor_cells For semi-supervised learning. Maximum number of anchor
#'   cells to use for each cell type.
#' @param min_anchor_cosine For semi-supervised learning. Cells must have at
#'   least this much cosine similarity to a fixed profile to be used as an
#'   anchor.
#' @param min_anchor_llr For semi-supervised learning. Cells must have
#'   (log-likelihood ratio / totalcounts) above this threshold to be used as an
#'   anchor
#' @param insufficient_anchors_thresh Cell types that end up with fewer than
#'   this many anchors will be discarded.
#' @param refinement Logical, flag for further anchor refinement via UMAP projection (default = FALSE)
#' @param blacklist vector of genes to be excluded for cell typing (default = NULL)
#' @param rescale Logical, flag for platform effect correction (default = FALSE).
#' @param refit Logical, flag for fitting reference profiles to anchors, run after rescale if rescale = TRUE (default = TRUE)
#' @return a list 
#' \describe{
#'     \item{updated_profiles}{a genes * cell types matrix for final updated reference profiles}
#'     \item{blacklist}{a vector of genes excluded from the final updated reference profiles}
#'     \item{anchors}{a named vector for final anchors used for reference profile update}
#'     \item{rescale_res}{a list of 5 elements, `rescaled_profiles`, `platformEff_statsDF`, `anchors`, `blacklist` and `lostgenes`, for platform effect correction outputs, return when rescale = TRUE}
#'     \item{refit_res}{a list of 2 elements, `refitted_profiles` and `anchors`, for anchor-based profile refitting outputs, return when refit = TRUE}
#' }
#' @export
#' @examples
#' data("mini_nsclc")
#' data("ioprofiles")
#' counts <- mini_nsclc$counts
#' ## estimate per-cell bg as a fraction of total counts:
#' negmean.per.totcount <- mean(rowMeans(mini_nsclc$neg)) / mean(rowSums(counts))
#' per.cell.bg <- rowSums(counts) * negmean.per.totcount
#' astats <- get_anchor_stats(counts = mini_nsclc$counts, 
#'                            assay_type="RNA", 
#'                            neg = Matrix::rowMeans(mini_nsclc$neg),
#'                            profiles = ioprofiles,
#'                            sds=NULL)
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
#'                                     insufficient_anchors_thresh = 5,
#'                                     assay_type="RNA")
#'
#' # The next step is to use the anchors to update the reference profiles:
#'
#' updateReferenceProfiles(reference_profiles = ioprofiles,
#'                         reference_sds = NULL,
#'                         counts = mini_nsclc$counts, 
#'                         neg = mini_nsclc$neg, 
#'                         assay_type = "rna", 
#'                         bg = per.cell.bg,
#'                         anchors = anchors) 

updateReferenceProfiles <-
  function(reference_profiles,
           reference_sds,
           counts,
           neg,
           assay_type = c("rna", "protein"), 
           bg = NULL,
           nb_size = 10,
           anchors = NULL,
           n_anchor_cells = 2000,
           min_anchor_cosine = 0.3,
           min_anchor_llr = 0.01,
           insufficient_anchors_thresh = 20,
           refinement = FALSE, 
           blacklist = NULL,
           rescale = FALSE, 
           refit = TRUE) {
    assay_type <- match.arg(tolower(assay_type), c("rna", "protein"))
    
    if(!any(rescale, refit)){
      stop("At least one of `rescale` or `refit` must be TRUE to update the reference profiles.")
    }
    
    # align genes:
    sharedgenes <- intersect(colnames(counts), rownames(reference_profiles))
    if(!is.null(blacklist)){
      sharedgenes <- setdiff(sharedgenes, blacklist)
    }
    
    ## step 1: if no initial anchors are provided, select them automatically:
    if (is.null(anchors)) {
      message("Automatically selecting anchor cells with the best fits to fixed profiles.")
      
      anchors <- find_anchor_cells(counts = counts[, sharedgenes], 
                                   neg = neg, 
                                   bg = bg, 
                                   profiles = reference_profiles[sharedgenes, ], 
                                   sds = reference_sds[sharedgenes, ], 
                                   size = nb_size, 
                                   assay_type = assay_type,
                                   n_cells = n_anchor_cells, 
                                   min_cosine = min_anchor_cosine, 
                                   min_scaled_llr = min_anchor_llr,
                                   insufficient_anchors_thresh = insufficient_anchors_thresh, 
                                   align_genes = FALSE, 
                                   refinement = refinement) 
    } else {
      # test anchors are valid:
      if(is.null(names(anchors))){
        if (length(anchors) != nrow(counts)) {
          stop("anchors must have length equal to the number of cells (row) in counts")
        }
        names(anchors) <- rownames(counts)
      } else {
        anchors <- anchors[rownames(counts)]
        
        if (all(is.na(anchors))) {
          stop("The provided `anchors` have no shared cell_id in names as the rownames in `counts`.")
        } 
      }
      
      
      # refine existing anchors 
      if (refinement){
        anchors <- refineAnchors(counts = counts[, sharedgenes], 
                                 neg = neg, 
                                 bg = bg, 
                                 align_genes = FALSE,
                                 profiles = reference_profiles[sharedgenes, ],  
                                 anchor_candidates = anchors, 
                                 nn_cells = n_anchor_cells,
                                 insufficient_anchors_thresh = insufficient_anchors_thresh)
      }
    }
    
    if (is.null(anchors))  {
      stop("No anchors were selected. The algorithm can't run under these conditions. 
         Solutions include: 1. make anchor selection more generous, try `refinement = FALSE`. 2. select anchors by hand.")
    }
    
    outs <- list()
    
    
    ## step 2: rescale profiles for platform effect based on high confidence of anchors
    if(rescale){
      message("Rescale reference profiles for platform effect.")
      outs[['rescale_res']] <- estimatePlatformEffects(counts = counts[, sharedgenes], 
                                                       neg = neg, 
                                                       bg = bg, 
                                                       assay_type = assay_type,
                                                       anchors = anchors,
                                                       profiles = reference_profiles[sharedgenes, ],
                                                       sds = reference_sds[sharedgenes, ])
      
      # add outliers from platform effect estimation, but exclude lostgenes
      blacklist <- unique(c(blacklist, outs[['rescale_res']][['blacklist']]))
      if(!is.null(blacklist)){
        sharedgenes <- setdiff(sharedgenes, blacklist)
      }
      
      if(!is.null(outs[['rescale_res']][['lostgenes']])){
        # add lostgenes with beta =1  back to updated profiles 
        outs[['updated_profiles']] <- rbind(outs[['rescale_res']][['rescaled_profiles']], 
                                            reference_profiles[outs[['rescale_res']][['lostgenes']], ]) 
        outs[['updated_sds']] <- rbind(outs[['rescale_res']][['rescaled_sds']], 
                                            reference_sds[outs[['rescale_res']][['lostgenes']], ]) 
      }else{
        outs[['updated_profiles']] <- outs[['rescale_res']][['rescaled_profiles']]
        outs[['updated_sds']] <- outs[['rescale_res']][['rescaled_sds']]
        
      }
    }
    
    
    ## step 3: second around of anchor selection based on corrected profiles
    if (rescale && refit){
      message("Second round of anchor selection given the rescaled reference profiles (lostgenes included with Beta =1): ")
      anchors_second <- find_anchor_cells(counts = counts[, sharedgenes], 
                                          neg = neg, 
                                          bg = bg, 
                                          profiles = outs[['updated_profiles']][sharedgenes, ], 
                                          sds = outs[['updated_sds']][sharedgenes, ], 
                                          size = nb_size, 
                                          n_cells = n_anchor_cells, 
                                          min_cosine = min_anchor_cosine, 
                                          min_scaled_llr = min_anchor_llr,
                                          insufficient_anchors_thresh = insufficient_anchors_thresh, 
                                          assay_type=assay_type,
                                          align_genes = FALSE, 
                                          refinement = refinement) 
      if (is.null(anchors_second)){
        message("No anchors were selected in the second around. Fall back to use first round of anchors for refitting. 
                Consider to 1. make anchor selection more generous, try `refinement = FALSE`. 2. select anchors by hand. 3. do `rescale` and `refit` in separate step. ")
      } else {
        # combine both rounds of anchors, which are in same cell_id orders as counts 
        anchors_second <- anchors_second[names(anchors)]
        combined_anchors <- sapply(
          seq_along(anchors), 
          function(idx){
            ct <- setdiff(unique(c(anchors[idx], anchors_second[idx])), NA)
            if(length(ct)!=1){
              # not to use if conflicts btw 2 rounds of anchor selection
              return(NA)
            }else {
              return(ct)
            }
          })
        names(combined_anchors) <- names(anchors)
        anchors <- combined_anchors
      }
    }
    
    
    ## step 4: refit the reference profiles using the second around of anchors 
    if(refit){
      # refit original reference profiles given the anchors, will include lostgenes from platform effect estimation  
      refitted_profiless <- updateProfilesFromAnchors(counts = counts[, sharedgenes],  
                                                     neg = neg, 
                                                     bg = bg,
                                                     anchors = anchors,
                                                     assay_type=assay_type)
      outs[['refit_res']] <- list(refitted_profiles = refitted_profiless$updated_profiles, 
                                  anchors = anchors)
      outs[['updated_profiles']] <- refitted_profiless$updated_profiles
      outs[['updated_sds']] <- refitted_profiless$updated_sds
    }
    
    outs[['blacklist']] <- blacklist
    outs[['anchors']] <- anchors
    
    return(outs)
  }



#' Use anchor cells to update reference profiles, simply by taking the mean
#' profile of the anchors.
#'
#' Uses anchor cells to estimate platform effects / scaling factors to be
#' applied to the genes/rows of the reference profile matrix. Then uses Bayesian
#' math to update the individual elements on X.
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell. Can be provided
#' @param bg Expected background
#' @param anchors Vector of anchor assignments
#' @param assay_type Assay type of RNA, protein (default = "rna")
#' 
#' @return \enumerate{ 
#' \item updated_profiles: A mean profiles matrix with the rows rescaled
#' according to platform effects and individual elements updated further 
#' \item updated_sds: A mean profiles matrix with the rows rescaled
#' according to platform effects and individual elements updated further}
#' @export
updateProfilesFromAnchors <-
  function(counts,
           neg,
           bg = NULL, 
           assay_type = c("rna", "protein"),
           anchors) {
    assay_type <- match.arg(tolower(assay_type), c("rna", "protein"))
    
    bg <- estimateBackground(counts, neg, bg)
    use <- !is.na(anchors)
    updated_profiles_info <- Estep(counts = counts[use, ],
                          clust = anchors[use],
                          neg = bg[use], 
                          assay_type=assay_type)
    updated_profiles <- updated_profiles_info$profiles
    updated_sds <- updated_profiles_info$sds
    
    return(list(updated_profiles=updated_profiles, 
                updated_sds=updated_sds))
  }


#' Platform effect adjustment on reference profiles based on the expression profiles of anchors 
#' 
#' Calculates gene-wise scaling factor between reference profiles and the observed profiles of the provided anchors.
#' @param counts Counts matrix, cells * genes.
#' @param neg Vector of mean negprobe counts per cell
#' @param assay_type Assay type of RNA, protein (default = "rna")
#' @param bg Expected background
#' @param anchors Vector giving "anchor" cell types, for use in semi-supervised
#'   clustering. Vector elements will be mainly NA's (for non-anchored cells)
#'   and cell type names for cells to be held constant throughout iterations.
#' @param profiles Matrix of reference profiles holding mean expression of genes x cell types. 
#'  Input linear-scale expression, with genes in rows and cell types in columns.
#' @param sds Matrix of reference profiles holding SDs expression of genes x cell types. 
#'  Input linear-scale expression, with genes in rows and cell types in columns. Only for assay_type of protein
#' @param blacklist vector of user-defined genes to be excluded for cell typing (default = NULL)
#' @return A list with five elements: 
#' \describe{
#'     \item{rescaled_profiles}{genes * cell types Matrix of rescaled reference profiles with platform effect corrected }
#'     \item{platformEff_statsDF}{a data.frame for statistics on platform effect estimation with genes in rows and columns for `Gene`, `Beta`, `beta_SE`.}
#'     \item{anchors}{a named vector of anchors used for platform effect estimation}
#'     \item{blacklist}{a vector of genes excluded from cell typing, including both outliers identified in platform effect estimation and the user-defined genes}
#'     \item{lostgenes}{a vector of genes excluded from platform effect estiamtion due to insufficient evidence}
#' }
#' @description The general workflow would be: (1) extract the anchor cells from input; 
#' (2) Run poisson regression with anchor cells; (3) Filter user defined genes(if any) 
#' and genes with extreme betas, outside [0.01, 100]; (4) Re-scale Profile with Beta estimates. 
#' @importFrom data.table data.table melt
#' @importFrom Matrix t rowSums
#' @export
estimatePlatformEffects <- 
  function(counts,
           neg, 
           assay_type = c("rna", "protein"),
           bg=NULL, 
           anchors, 
           profiles,
           sds=NULL,
           blacklist=NULL){
    assay_type <- match.arg(tolower(assay_type), c("rna", "protein"))
    
    #### Step1: clean up and prepare inputs, focus on anchors only
    bg <- estimateBackground(counts = counts, neg = neg, bg = bg)
    # non-NA anchors 
    anchors <- anchors[!is.na(anchors)]
    anchors <- anchors[(anchors %in% colnames(profiles)) & (names(anchors) %in% rownames(counts))]
    cts_to_check <- unique(anchors)
    
    if(length(anchors)<20 || length(cts_to_check)<3){
      stop(sprintf("Only %d anchor for %d cell types shared among `counts` and `profiles`, check if missing names for cell_id in `anchors`.", 
                   length(anchors), length(cts_to_check)))
    } else {
      message(sprintf("Estimate platform effect given %d anchor across %d cell types.", 
                      length(anchors), length(cts_to_check)))
    }
    
    # focus on shared genes 
    counts <- alignGenes(counts = counts[names(anchors), ], 
                             profiles = profiles)
    profiles <- profiles[colnames(counts), ]
    sds <- sds[colnames(counts), ]
    
    ## group shared genes based on their expression level in obs vs. ref
    # (1) ok ref & above-zero obs, evaluate in glm; 
    # (2) ok ref but near-zero obs, add to blacklist as outliers; 
    # (3) near-zero ref but high obs, add to blacklist as outliers;
    # (4) near-zero ref but low obs, add to lostgenes. 
    
    # dense array for net count, high memory consumption  
    if(identical(tolower(assay_type), "rna")){
      netCount_data <- pmax(counts - bg[names(anchors)], 0)
    }
    if(identical(tolower(assay_type), "protein")){
      netCount_data <- counts 
    }
    netAvg_perCT <- sapply(cts_to_check, function(ct){
      Matrix::colMeans(netCount_data[anchors == ct, , drop = F], na.rm = T)
    })
    gc()
    
    geneDF <- data.frame(
      Gene = colnames(counts), 
      Ref_maxAll = apply(profiles[, cts_to_check, drop = F], 1, max, na.rm = T),
      NetCount_maxAll = apply(netCount_data, 2, max, na.rm = T),
      NetCount_maxPerCT = apply(netAvg_perCT, 1, max, na.rm = T)
    )
    rm(netCount_data, netAvg_perCT)
    gc()
    
    # ok ref > 0.5* median(ref of all anchor cell types)
    geneDF[['ok_ref']] <- (geneDF[['Ref_maxAll']] > 0.5* median(profiles[, cts_to_check, drop = F], na.rm = T))
    # positive obs, max (net of all) >=1
    geneDF[['pos_net']] <- (geneDF[['NetCount_maxAll']] >= 1)
    # high obs, max(average per cell type) >=2
    geneDF[['high_net']] <- (geneDF[['NetCount_maxPerCT']] >= 2)
    
    # lostgenes: insufficient evidence for evaluation 
    lostgenes <- geneDF$Gene[!(geneDF$ok_ref | geneDF$high_net)]
    if(length(lostgenes)>0){
      message(paste0(length(lostgenes), " genes with low expression in both `counts` and `profiles` are excluded from platform estimation."))
    }else{
      lostgenes <- NULL
    }
    
    # blacklist with contradictory trend
    blacklist_addon <- geneDF$Gene[(geneDF$ok_ref & !geneDF$pos_net) | (!geneDF$ok_ref & geneDF$high_net)]
    if(length(blacklist_addon)>0){
      message(paste0(length(blacklist_addon), " genes with contradictory trend are appended to `blacklist`. "))
      blacklist <- unique(c(blacklist, blacklist_addon))
    }
    
    # genes for evaluation 
    use_genes <- geneDF$Gene[geneDF$ok_ref & geneDF$pos_net]
    counts <- counts[,use_genes, drop = F]
    profiles <- profiles[use_genes, ]
    sds <- sds[use_genes, ]
    if(ncol(counts)<1){
      stop("No shared genes with sufficient count for platform evaluation.")
    }
    
    
    # fast conversion to data.frame for glm
    count_vector <- data.table::melt(data.table::data.table(as.matrix(counts), keep.rownames = TRUE) , id.vars = c("rn"))
    colnames(count_vector) <- c("Cell_ID","Gene","Counts")
    
    # expanded cell types * genes matrix 
    reference_data <- t(as.matrix(profiles[, unname(anchors)]))
    reference_vector <- data.table::melt(data.table::data.table(reference_data, keep.rownames = TRUE) , id.vars = c("rn"))
    colnames(reference_vector) <- c("Cell_Type","Gene","Reference")
    
    ### Step2: Run poisson regression with link being Identity to estimate gene-level platform effect based on selected anchor cells
    query_DF <- as.data.frame(cbind(count_vector, reference_vector[, Gene:=NULL]))
    query_DF[['BG']] <- bg[count_vector$Cell_ID]
    # cell level scaling factor between obs vs. reference 
    query_DF[['Cell_SF']] <- Matrix::rowSums(counts)/Matrix::rowSums(reference_data)
    
    # split by gene
    query_DF <- split(query_DF, query_DF$Gene)
    
    # fastglm estimation with method = 3L to allow non-positive definite matrices 
    percentCores <- 0.25
    PlatformEstimator_fastGLM <- function(df){
      if(identical(tolower(assay_type), "rna")){
        GLM_Fit<- fastglm::fastglm(x = model.matrix(~Reference : Cell_SF -1, data = df), 
                                   y = df$Counts, offset = df$BG, 
                                   family= stats::poisson(link = "identity"),
                                   start=c(0), method = 3L)
      }
      if(identical(tolower(assay_type), "protein")){
        GLM_Fit<- fastglm::fastglm(x = model.matrix(~Reference : Cell_SF -1, data = df), 
                                   y = df$Counts, # offset = df$BG, 
                                   family= stats::gaussian(link = "identity"),
                                   start=c(0), method = 3L)
      }

      return(data.frame(Gene=df$Gene[1],
                        Beta=GLM_Fit$coefficients["Reference:Cell_SF"],
                        beta_SE=summary(GLM_Fit)$coefficients["Reference:Cell_SF","Std. Error"]))
    }
    PlatformEff <- parallel::mclapply(query_DF, PlatformEstimator_fastGLM, 
                                      mc.cores = numCores(percentCores =percentCores))
    
    # give warnings for estimation error
    warning_genes <- (sapply(PlatformEff, class) != "data.frame")
    if(sum(warning_genes)>0){
      message(paste0(sum(warning_genes), " genes with error in glm fitting are appended to `lostgenes`, try set option with smaller `mc.cores`."))
      lostgenes <- unique(c(lostgenes, names(warning_genes)[which(warning_genes)]))
      PlatformEff <- PlatformEff[!warning_genes]
    }
    PlatformEff <- as.data.frame(do.call(rbind, PlatformEff))
    rownames(PlatformEff) <- PlatformEff$Gene
    
    
    ### Step3: Filtering genes with extreme betas, outside [0.01, 100] and pre-specified by the user
    blacklist_addon <- as.vector(PlatformEff$Gene[PlatformEff$Beta < 0.01 | PlatformEff$Beta > 100])
    if(length(blacklist_addon)>0){
      message(paste0(length(blacklist_addon), " genes with extreme beta are appended to `blacklist`. "))
      blacklist <- unique(c(blacklist, blacklist_addon))
    }
    genes_to_keep <- setdiff(PlatformEff$Gene, blacklist)
    
    ### Step4: Rescale the raw reference profile
    rownames(PlatformEff) <- PlatformEff$Gene
    rescaled_profiles <- diag(PlatformEff[genes_to_keep,]$Beta) %*% profiles[genes_to_keep,]
    rownames(rescaled_profiles) <- genes_to_keep
    
    if(identical(tolower(assay_type), "protein")){
      rescaled_sds <- diag(PlatformEff[genes_to_keep,]$Beta) %*% sds[genes_to_keep,]
      rownames(rescaled_sds) <- genes_to_keep
    }else{
      rescaled_sds <- sds[genes_to_keep,]
    }

    
    
    return(list(rescaled_profiles = rescaled_profiles,
                rescaled_sds = rescaled_sds,
                platformEff_statsDF = PlatformEff,
                anchors = anchors,
                blacklist = blacklist,
                lostgenes = lostgenes))
    
  }

