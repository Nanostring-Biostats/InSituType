# load data ("raw" and "cellannot"):
data("ioprofiles")
data("iocolors")
data("mini_nsclc")

# test supervised cell typing:
sup <- insitutype(counts = mini_nsclc$counts,
                  neg = Matrix::rowMeans(mini_nsclc$neg),
                  bg = NULL,
                  init_clust = NULL, 
                  n_clusts = 0,
                  fixed_profiles = ioprofiles[,1:6],
                  nb_size = 10,
                  method = "EM",
                  n_starts = 2,
                  align_genes = TRUE,
                  n_benchmark_cells = 200,
                  n_phase1 = 500,
                  n_phase2 = 1000,
                  n_phase3 = 2000,
                  pct_drop = 1/10000, 
                  min_prob_increase = 0.05,
                  max_iters = 10)   

testthat::test_that("supervised cell typing produces correct outputs", {
  expect_true(all(is.element(c("clust", "probs"), names(sup))))
  expect_true(is.vector(sup$clust))
  expect_true(is.matrix(sup$probs))
})


# run unsupervised clustering with several random starts:
unsup <- insitutype(counts = mini_nsclc$counts,
                    neg = Matrix::rowMeans(mini_nsclc$neg),
                    bg = NULL,
                    init_clust = NULL, n_clusts = 2:5,
                    fixed_profiles = NULL,
                    nb_size = 10,
                    method = "EM",
                    n_starts = 2,
                    align_genes = TRUE,
                    sketchingdata = NULL,
                    n_benchmark_cells = 100,
                    n_phase1 = 50,
                    n_phase2 = 100,
                    n_phase3 = 200,
                    n_chooseclusternumber = 100,
                    pct_drop = 1/10000, 
                    min_prob_increase = 0.05,
                    max_iters = 10)   

testthat::test_that("unsupervised cell typing produces correct outputs", {
  expect_true(all(is.element(c("clust", "probs", "profiles"), names(unsup))))
  expect_true(is.vector(unsup$clust))
  expect_true(is.matrix(unsup$probs))
  expect_true(is.matrix(unsup$profiles))
})




# run unsupervised clustering with init_clust specified:
init_clust <- rep(c("name1", "xxx", "ooo"), each = nrow(mini_nsclc$counts) / 3)[1:nrow(mini_nsclc$counts)]
unsup <- insitutype(counts = mini_nsclc$counts,
                    neg = Matrix::rowMeans(mini_nsclc$neg),
                    bg = 0.03,
                    init_clust = init_clust, 
                    n_clusts = 6,
                    fixed_profiles = NULL,
                    nb_size = 10,
                    method = "EM",
                    n_starts = NULL,
                    align_genes = TRUE,
                    sketchingdata = NULL,
                    n_benchmark_cells = 100,
                    n_phase1 = 50,
                    n_phase2 = 100,
                    n_phase3 = 200,
                    pct_drop = 1/10000, 
                    min_prob_increase = 0.05,
                    max_iters = 10)   


testthat::test_that("unsupervised cell typing using init_clust produces correct outputs", {
  expect_true(all(is.element(c("clust", "probs", "profiles"), names(unsup))))
  expect_true(is.vector(unsup$clust))
  expect_true(is.matrix(unsup$probs))
  expect_true(is.matrix(unsup$profiles))
  expect_true(all(sort(unique(unsup$clust)) == sort(unique(init_clust))))
})

# plot clusters:
if (FALSE) {
  scols = c(iocolors, 
            '#8DD3C7','#BEBADA','#FB8072','#80B1D3','#FDB462','#B3DE69','#FCCDE5','#D9D9D9','#BC80BD','#CCEBC5','#FFED6F','#66C2A5','#FC8D62','#8DA0CB','#E78AC3','#A6D854','#FFD92F','#E5C494','#B3B3B3')
  scols <- scols[1:length(unique(unsup$clust))]
  names(scols) = unique(unsup$clust)
  plot(mini_nsclc$x, mini_nsclc$y, pch = 16, col = scols[unsup$clust])
}

# view immune oncology cell profiles (in ptolemy package data):
semi <- insitutype(counts = mini_nsclc$counts,
                   neg = Matrix::rowMeans(mini_nsclc$neg),
                   bg = NULL,
                   init_clust = NULL, n_clusts = 2,
                   fixed_profiles = ioprofiles[, 1:3],
                   nb_size = 10,
                   method = "EM",
                   n_starts = 2,
                   align_genes = TRUE,
                   n_benchmark_cells = 200,
                   n_phase1 = 500,
                   n_phase2 = 1000,
                   n_phase3 = 2000,
                   n_chooseclusternumber = 100,
                   pct_drop = 1/10000, 
                   min_prob_increase = 0.05,
                   max_iters = 10)   

# run semisupervised clustering with init_clust specified:
init_clust <- rep(letters[1:3], each = nrow(mini_nsclc$counts) / 3)[1:nrow(mini_nsclc$counts)]
semi <- insitutype(counts = mini_nsclc$counts,
                   neg = Matrix::rowMeans(mini_nsclc$neg),
                   bg = NULL,
                   init_clust = init_clust, n_clusts = 3,
                   fixed_profiles = ioprofiles[, 1:3],
                   nb_size = 10,
                   method = "EM",
                   n_starts = 2,
                   align_genes = TRUE,
                   n_benchmark_cells = 200,
                   n_phase1 = 500,
                   n_phase2 = 1000,
                   n_phase3 = 2000,
                   pct_drop = 1/10000, 
                   min_prob_increase = 0.05,
                   max_iters = 10)   