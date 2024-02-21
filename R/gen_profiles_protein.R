#' Generate the mean reference profile and its SD reference profile based on the data itself
#' This function is based on signature matrix included in CELESTA package 
#' First, we rebuild a nested cell typing lists based on the 2-D signature matrix
#' Second, we identify anchor cells ranked by their expression level for each cell type's protein marker
#' Third, we estimate averaged expression level and SDs for proteins and cell types using the anchors
#'
#' @param exp.mat a matrix of raw protein expression data. cells are in rows and proteins are in columns
#' @param sig_mat a signature matrix of cell types. cell types x protein markers 
#' @param cutoff a cutoff of quantile. e.g) cutoff=0.9 means that top 90 percentiles of cells are called anchors for the protein expression
#' @param min.num.cells a minimum number of cells each cell type to estimate its mean or SDs. default value is 30.
#' @param keep_marker_proteins whether just marker proteins from the signature matrix is kept. default value is FALSE, which returns all proteins included in the data
#' 
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr summarise_all group_by filter
#' @return A list, with the following elements:
#' \enumerate{
#' \item mean.ref.profile: a matrix of cluster-specific expression profiles. proteins x cell types
#' \item SDs.ref.profile: a matrix of standard deviation profiles of pre-defined clusters. proteins x cell types
#' \item anchors: a vector giving "anchor" cell types. Vector elements will be mainly NA's (for non-anchored cells)
#' }
#' @name gen_profiles_protein_expression
#' @examples 
#' data("tonsil_protein")
#' data("human_signature")
#' data("mouse_signature")
#' references <- gen_profiles_protein_expression(
#'  exp.mat=tonsil_protein$counts,
#'  sig_mat=NULL)
gen_profiles_protein_expression <- function(exp.mat, sig_mat=NULL, cutoff=0.9, min.num.cells=30, keep_marker_proteins=FALSE){

  if(is.null(sig_mat)){

    ## call the human's signature matrix
    sig_mat = InSituType::human_signature
    
    ## If the panel is for mouse, we call the mouse's signature matrix
    if(length(intersect(names(sig_mat), names(exp.mat)) == 0)){
      sig_mat = InSituType::mouse_signature
    }
  }

  markerProteins <- intersect(colnames(sig_mat), colnames(exp.mat))
  ## Split Lineage levels into columns
  sig_mat[is.na(sig_mat)] <- 0
  sig_mat$level1 <- lapply(strsplit(sig_mat$Lineage_level, "_"), function(x){x[1]}) %>% unlist()
  sig_mat$level2 <- lapply(strsplit(sig_mat$Lineage_level, "_"), function(x){x[2]}) %>% unlist()
  sig_mat$level3 <- lapply(strsplit(sig_mat$Lineage_level, "_"), function(x){x[3]}) %>% unlist()
  
  markerProtein_celltype_level <- vector("list", length=max(sig_mat$level1))
  for (i in 1:max(sig_mat$level1)){
    
    if(i ==1){
      markerProtein_celltype_level[[i]] <- data.frame(celltype = sig_mat[sig_mat$level1==i,]$celltype, 
                                                      marker_protein=apply(sig_mat[sig_mat$level1==i,], 1, function(x){colnames(sig_mat[sig_mat$level1==i,])[which(x==1)[1]]}),
                                                      upper_celltype = "Parent")
    }else{
      markerProtein_celltype_level[[i]] <- data.frame(celltype = sig_mat[sig_mat$level1==i,]$celltype, 
                                                      marker_protein=apply(sig_mat[sig_mat$level1==i,], 1, function(x){colnames(sig_mat[sig_mat$level1==i,])[which(x==1)[1]]}),
                                                      upper_celltype = sig_mat$celltype[which(sig_mat$level3==unique(sig_mat[sig_mat$level1==i,]$level2))])
    }
  }
  
  dat_mat_level <- vector("list", length=max(sig_mat$level1))
  for (i in 1:length(markerProtein_celltype_level)){
    if(i ==1){
      dat_mat_level[[i]] <- lapply(markerProtein_celltype_level[[i]]$marker_protein, function(x){
        if(max(exp.mat)<=1){
          cutoff <- 0.9
        }else{
          cutoff <- quantile(exp.mat[, x], prob=0.9)
        }
        rownames(exp.mat)[which(exp.mat[, x] > cutoff)]
      })
      names(dat_mat_level[[i]]) <- markerProtein_celltype_level[[i]]$celltype
    }else{
      
      dat_mat_level[[i]] <- vector("list", nrow(markerProtein_celltype_level[[i]]))
      names(dat_mat_level[[i]]) <- markerProtein_celltype_level[[i]]$celltype
      
      for(j in 1:length(markerProtein_celltype_level[[i]]$celltype)){
        
        if(!is.na(markerProtein_celltype_level[[i]][j,]$marker_protein)){
          
          ## Identify the upper level's cell type and where it is located in the signature matrix' lineage level
          for(k in 1:(i-1)){
            tempDD <- markerProtein_celltype_level[[k]] %>% filter(celltype==markerProtein_celltype_level[[i]]$upper_celltype[1])
            
            if(nrow(tempDD)==1){
              tempMar=tempDD
              idx_k=k
            }else{
              paste("pass")
            }
          }
          
          tempD <- exp.mat[rownames(exp.mat) %in% dat_mat_level[[idx_k]][[tempMar$celltype]], ]
          
          if(max(exp.mat)<=1){
            cutoff <- 0.9
          }else{
            cutoff <- quantile(tempD[, markerProtein_celltype_level[[i]][j,]$marker_protein], prob=0.9)
          }
          
          tempID <- rownames(tempD)[which(tempD[, markerProtein_celltype_level[[i]][j,]$marker_protein] > cutoff)]
          
          dat_mat_level[[idx_k]][[tempMar$celltype]] <- setdiff(dat_mat_level[[idx_k]][[tempMar$celltype]], tempID)       
          dat_mat_level[[i]][[markerProtein_celltype_level[[i]]$celltype[j]]] <- tempID
        }else{
          break
          
        }
      }
    }
  }
  
  markerProtein_celltype_all <- do.call("rbind", markerProtein_celltype_level)
  marker_id_cell_type <- do.call(c, dat_mat_level)
  marker_id_cell_type_insitu <- marker_id_cell_type[lapply(marker_id_cell_type, length)!=0]
  marker_id_cell_type_insitu_df <- lapply(1:length(marker_id_cell_type_insitu), function(x){data.frame(cell_ID=marker_id_cell_type_insitu[[x]], 
                                                                                                    celltype=rep(names(marker_id_cell_type_insitu[x]), length(marker_id_cell_type_insitu[[x]])))})
  names(marker_id_cell_type_insitu_df) <- names(marker_id_cell_type_insitu)
  anchors <- do.call("rbind", marker_id_cell_type_insitu_df) %>% as.data.frame()
  anchors_duplicate <- anchors[which(duplicated(anchors$cell_ID)==TRUE),]$cell_ID
  
  marker_id_cell_type_unique <- lapply(marker_id_cell_type_insitu_df, 
                                       function(x) {
                                         tempV <- setdiff(x$cell_ID, anchors_duplicate)
                                         if(length(tempV) > 20){
                                           names(tempV) <- x[x$cell_ID %in% tempV,]$celltype
                                           tempV <- tempV
                                         }else{
                                           tempV <- NULL
                                         }
                                         return(tempV)})
  
  # marker_id_cell_type_unique <- Filter(Negate(is.null), marker_id_cell_type_unique)
  
  anchors <- anchors[which(duplicated(anchors$cell_ID)==FALSE),] 
  anchors <- anchors %>% filter(celltype %in% names(marker_id_cell_type_unique))
  
  anchors <- rbind(anchors, data.frame(cell_ID = setdiff(rownames(exp.mat), anchors$cell_ID), celltype=NA))
  rownames(anchors) <- anchors$cell_ID
  anchors$cell_ID <- NULL
  anchors <- t(anchors)[1,]
  
  ############################ Estimate averaged protein expression each cell type with its anchor cells ######################################
  protein_exp_means_list <- lapply(marker_id_cell_type_unique, function(x){
    
    mean.exp <- exp.mat[rownames(exp.mat) %in% x, ] %>% colMeans()
    
  })
  
  mean.ref.profile <- do.call("rbind", protein_exp_means_list) %>% t() %>% as.data.frame()
  
  protein_exp_SDs_list <- lapply(marker_id_cell_type_unique, function(x){
    apply(exp.mat[rownames(exp.mat) %in% x, ], 2, sd )
  })
  names(protein_exp_SDs_list) <- names(marker_id_cell_type_unique)
  SDs.ref.profile <- do.call("rbind", protein_exp_SDs_list) %>% t() %>% as.data.frame()
  
  if(keep_marker_proteins){
    mean.ref.profile <- mean.ref.profile[markerProteins, ]
    SDs.ref.profile <- SDs.ref.profile[markerProteins, ]
  }
  out <- list(mean.ref.profile=mean.ref.profile, SDs.ref.profile=SDs.ref.profile, anchors=anchors[rownames(exp.mat)])
  return(out)
}


#' Generate the mean reference profile and its SD reference profile from an annotation file
#' This function is only for protein data set with known anchor cells and their cell types
#'
#' @param exp.mat a matrix of raw protein expression data. cells are in rows and proteins are in columns
#' @param anno a data frame or matrix of cell types for anchor cells or manually annotated cell typing information for some cells. Should include cell_ID and celltype at least. 
#' 
#' @return A list, with the following elements:
#' \enumerate{
#' \item mean.ref.profile: a matrix of cluster-specific expression profiles. proteins * cell types
#' \item SDs.ref.profile: a matrix of standard deviation profiles of pre-defined clusters. proteins * cell types
#' \item anchors: a vector giving "anchor" cell types. Vector elements will be mainly NA's (for non-anchored cells)
#' }

gen_profiles_protein_annotation <- function(exp.mat, anno) {
  
  anno_ref_mat <- merge(exp.mat %>% as.data.frame() %>% rownames_to_column(var="cell_ID"), anno %>% dplyr::select(c(cell_ID, cellType)), by="cell_ID") %>% column_to_rownames(var="cell_ID")
  
  mean.ref.profile <- anno_ref_mat %>% group_by(cellType) %>% summarise_all(mean) %>% column_to_rownames(var="cellType") %>% t()
  SDs.ref.profile <- anno_ref_mat %>% group_by(cellType) %>% summarise_all(sd) %>% column_to_rownames(var="cellType") %>% t()
  
  ## Set NAs for non-anchor cells' cell types
  anchors <- rbind(anno %>% dplyr::select(c(cell_ID, cellType)), data.frame(cell_ID = setdiff(rownames(exp.mat), rownames(anno_ref_mat)), cellType=NA)) 
  rownames(anchors) <- anchors$cell_ID
  anchors$cell_ID <- NULL
  anchors <- anchors %>% t()
  anchors <- anchors[1,]
  
  out <- list(mean.ref.profile=mean.ref.profile,
              SDs.ref.profile=SDs.ref.profile,
              anchors=anchors[rownames(exp.mat)])
  return(out)
}
