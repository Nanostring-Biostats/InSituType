library(testthat)

data("mini_tma")

counts <- t(as.matrix(mini_tma@expression$rna$raw))
s <- rowSums(counts)
neg <- Matrix::colMeans(mini_tma@expression$neg$raw)

# run a very fast example:
res <- chooseClusterNumber(counts = counts, s = s, neg = neg, 
                           n_clusts = 2:4, 
                           n_iters = 4, subset_size = 500)
  

# version with fixed profiles:
res2 <- chooseClusterNumber(counts = counts, s = s, neg = neg, 
                           n_clusts = 2:4, 
                           n_iters = 4, subset_size = 500,
                           fixed_profiles = ioprofiles)

# what to test:
# - a chosen cluster number is output, and it's from the N with the smallest AIC
