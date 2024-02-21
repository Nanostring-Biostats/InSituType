
#' Fast approximate PCA
#' 
#' Run fast approximate PCA by deriving components from a subset then applying to all observations
#' This is implemented to avoid version issues with irlba, which actually runs faster than this. 
#' @param mat Regular or sparse Matrix, observations in rows, variables in columns
#' @param ncomp Number of components to estimate
#' @param nsub Subsample size
#' @return A matrix of n x ncomp approximate PC values
fastApproxPCA <- function(mat, ncomp, nsub = 5000) {
  
  # define subset:
  if (nsub < nrow(mat)) {
    # deterministic but random-like subset: just sample at regular intervals
    n <- nrow(mat)
    nsub <- min(nrow(mat), nsub)
    sub <- setdiff(unique(((1:nsub) * round(n / nsub )) %% n), 0)
  } else {
    sub <- TRUE
  }
  
  # PCA on subset:
  pc <- prcomp(mat[sub, , drop = FALSE])$rotation[, 1:ncomp]
  
  # apply to whole matrix:
  return(mat %*% pc)
}

# confirming speed is OK:
if (FALSE) {
  mat = c()
  for (i in 1:20){
    mat = rbind(mat, mini_nsclc$counts)
  }
  
  tic()
  r1 = fastApproxPCA(mat, ncomp = 20, nsub = 5000) 
  toc()
  
  tic()
  r2 = irlba::prcomp_irlba(mat, n = 20) 
  toc()
  
  diag(cor(r1, r2$x))
  
}
