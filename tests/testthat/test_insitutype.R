set.seed(0)

# load data ("raw" and "cellannot"):
data("ioprofiles")
data("iocolors")
data("mini_nsclc")

if (FALSE) {
  library(testthat)
  library(MLEcell)
  nbclust = MLEcell:::nbclust
  Mstep=MLEcell:::Mstep
  Estep=MLEcell:::Estep
  runinsitutype=MLEcell:::runinsitutype
}

# for line-by-line debugging:
if (FALSE) {
  counts = mini_nsclc$counts;neg = Matrix::rowMeans(mini_nsclc$neg);bg = NULL;anchors = NULL
  init_clust = NULL; n_clusts = 2:3;fixed_profiles = ioprofiles[, 1:3];
  nb_size = 10;n_starts = 2;align_genes = TRUE;n_benchmark_cells = 200;n_phase1 = 300;n_phase2 = 500;n_phase3 = 1000;
  n_chooseclusternumber = 300;pct_drop = 1/5000;min_prob_increase = 0.05;max_iters = 4;n_anchor_cells = 20
  min_anchor_cosine = 0.3; min_anchor_llr = 0.1;sketchingdata=NULL;anchor_replacement_thresh=0.75;insufficient_anchors_thresh = 1
  sketchingdata = NULL; anchor_replacement_thresh = 5
}

# test supervised cell typing using direct loglik calcs:
sup <- insitutypeML(counts = mini_nsclc$counts,
                      neg = Matrix::rowMeans(mini_nsclc$neg),
                      bg = NULL,
                      fixed_profiles = ioprofiles[,1:6],
                      nb_size = 10, 
                      align_genes = TRUE) 
  
testthat::test_that("supervised cell typing produces correct outputs", {
  expect_true(all(is.element(c("clust", "prob", "logliks"), names(sup))))
  expect_true(is.vector(sup$clust))
  expect_true(is.vector(sup$prob))
  expect_true(is.matrix(sup$logliks))
  
})

# test semi-supervised with 0 new clusts:
semi <- insitutype(counts = mini_nsclc$counts,
                  neg = Matrix::rowMeans(mini_nsclc$neg),
                  tissue = sample(letters[1:2], nrow(mini_nsclc$counts), replace = TRUE),
                  bg = NULL,
                  init_clust = NULL, 
                  n_clusts = 0,
                  anchors = NULL,
                  fixed_profiles = ioprofiles[,1:6],
                  nb_size = 10,
                  n_starts = 2,
                  align_genes = TRUE,
                  n_benchmark_cells = 200,
                  n_phase1 = 500,
                  n_phase2 = 1000,
                  n_phase3 = 2000,
                  pct_drop = 1/10000, 
                  min_prob_increase = 0.05,
                  max_iters = 4,
                  n_anchor_cells = 20, min_anchor_cosine = 0.3, min_anchor_llr = 0.01, insufficient_anchors_thresh = 2)   

testthat::test_that("semiservised cell typing with n_clusts = 0 produces correct outputs", {
  expect_true(all(is.element(c("clust", "prob", "logliks"), names(semi))))
  expect_true(is.vector(semi$clust))
  expect_true(is.vector(semi$prob))
  expect_true(is.matrix(semi$logliks))
})


# run unsupervised clustering with several random starts:
unsup <- insitutype(counts = mini_nsclc$counts,
                    neg = Matrix::rowMeans(mini_nsclc$neg),
                    bg = NULL,
                    tissue =  sample(letters[1:2], nrow(mini_nsclc$counts), replace = TRUE),
                    init_clust = NULL, n_clusts = 2:5,
                    fixed_profiles = NULL,
                    anchors = NULL,
                    nb_size = 10,
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
                    max_iters = 2)   

testthat::test_that("unsupervised cell typing produces correct outputs", {
  expect_true(all(is.element(c("clust", "prob", "logliks", "profiles"), names(unsup))))
  expect_true(is.vector(unsup$clust))
  expect_true(is.vector(unsup$prob))
  expect_true(is.matrix(unsup$logliks))
  expect_true(is.matrix(unsup$profiles))
})




# run unsupervised clustering with init_clust specified:
init_clust <- rep(c("name1", "xxx", "ooo"), each = nrow(mini_nsclc$counts) / 3)[1:nrow(mini_nsclc$counts)]
unsup <- insitutype(counts = mini_nsclc$counts,
                    neg = Matrix::rowMeans(mini_nsclc$neg),
                    bg = 0.03,
                    tissue = NULL,
                    init_clust = init_clust, 
                    n_clusts = 6,
                    fixed_profiles = NULL,
                    nb_size = 10,
                    n_starts = NULL,
                    align_genes = TRUE,
                    sketchingdata = NULL,
                    n_benchmark_cells = 100,
                    n_phase1 = 50,
                    n_phase2 = 100,
                    n_phase3 = 200,
                    pct_drop = 1/10000, 
                    min_prob_increase = 0.05,
                    max_iters = 4)   


testthat::test_that("unsupervised cell typing using init_clust produces correct outputs", {
  expect_true(all(is.element(c("clust", "prob", "logliks", "profiles"), names(unsup))))
  expect_true(is.vector(unsup$clust))
  expect_true(is.vector(unsup$prob))
  expect_true(is.matrix(unsup$logliks))
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

# semi-supervised using the immune oncology cell profiles (in ptolemy package data):
semi <- insitutype(counts = mini_nsclc$counts,
                   neg = Matrix::rowMeans(mini_nsclc$neg),
                   bg = NULL,
                   tissue =  sample(letters[1:2], nrow(mini_nsclc$counts), replace = TRUE),
                   anchors = NULL,
                   init_clust = NULL, n_clusts = 2,
                   fixed_profiles = ioprofiles[, 1:3],
                   nb_size = 10,
                   n_starts = 2,
                   align_genes = TRUE,
                   n_benchmark_cells = 200,
                   n_phase1 = 500,
                   n_phase2 = 1000,
                   n_phase3 = 2000,
                   n_chooseclusternumber = 100,
                   pct_drop = 1/5000, 
                   min_prob_increase = 0.05,
                   max_iters = 4,
                   n_anchor_cells = 20, 
                   min_anchor_cosine = 0.3, 
                   min_anchor_llr = 0.01,
                   sketchingdata = NULL, insufficient_anchors_thresh = 2)   


testthat::test_that("unsupervised cell typing using init_clust produces correct outputs", {
  expect_true(all(is.element(c("clust", "prob", "logliks", "profiles"), names(semi))))
  expect_true(is.vector(semi$clust))
  expect_true(is.vector(semi$prob))
  expect_true(is.matrix(semi$logliks))
  expect_true(is.matrix(semi$profiles))
})


# test merge cells with multi-sample clustering:
merge2 <- mergeCells(
  merges =  c("B-cell" = "lymphoid", "a_cl1" = "cancer", "b_cl2" = "cancer"), 
  to_delete = c("endothelial", "b_cl1"), 
  logliks = semi$logliks)
testthat::test_that("mergecells works when merges and deletions are asked for", {
  expect_true(all(is.element(colnames(merge2$logliks), c("lymphoid", "cancer", "fibroblast", "a_cl2"))))
  expect_true(all(is.na(merge2$logliks[semi$logliks[, "endothelial"] == 1, 1])))
  expect_equal(names(merge2$clust), names(semi$clust))
})





# semi supervised with choosing cluster number
semi <- insitutype(counts = mini_nsclc$counts,
                   neg = Matrix::rowMeans(mini_nsclc$neg),
                   bg = NULL,
                   tissue = NULL,
                   anchors = NULL,
                   init_clust = NULL, n_clusts = 2:3,
                   fixed_profiles = ioprofiles[, 1:3],
                   nb_size = 10,
                   n_starts = 2,
                   align_genes = TRUE,
                   n_benchmark_cells = 200,
                   n_phase1 = 300,
                   n_phase2 = 500,
                   n_phase3 = 1000,
                   n_chooseclusternumber = 300,
                   pct_drop = 1/5000, 
                   min_prob_increase = 0.05,
                   max_iters = 4,
                   n_anchor_cells = 20, min_anchor_cosine = 0.3, min_anchor_llr = 0.01, insufficient_anchors_thresh = 2)   


testthat::test_that("unsupervised cell typing using init_clust produces correct outputs", {
  expect_true(all(is.element(c("clust", "prob", "logliks", "profiles"), names(semi))))
  expect_true(is.vector(semi$clust))
  expect_true(is.vector(semi$prob))
  expect_true(is.matrix(semi$logliks))
  expect_true(is.matrix(semi$profiles))
})


# test chooseclusternumber:
res <- chooseClusterNumber(counts = mini_nsclc$counts,
                           neg = Matrix::rowMeans(mini_nsclc$neg),
                           bg = NULL,
                           anchors = NULL,
                           init_clust = NULL, n_clusts = 2:3,
                           fixed_profiles = ioprofiles[, 1:3],
                           max_iters = 4, 
                           subset_size = 1000, 
                           align_genes = TRUE, 
                           plotresults = FALSE, 
                           nb_size = 10, 
                           pct_drop = 0.005, 
                           min_prob_increase = 0.05) 

testthat::test_that("chooseClusterNumber produces correct outputs", {
  expect_true(all(is.element(c("best_clust_number", "n_clusts", "loglik", "aic", "bic"), names(res))))
  expect_true(length(res$best_clust_number) == 1)
  expect_true(length(res$n_clusts) == 2)
  expect_true(length(res$loglik) == 2)
  expect_true(length(res$aic) == 2)
  expect_true(length(res$bic) == 2)
})

# test anchor stats function:
astats <- get_anchor_stats(counts = mini_nsclc$counts,
                           neg = Matrix::rowMeans(mini_nsclc$neg),
                           bg = NULL, 
                           align_genes = TRUE,
                           profiles = ioprofiles[, 1:3], 
                           size = 10, 
                           min_cosine = 0.3) 

testthat::test_that("get_anchor_stats produces correct outputs", {
  expect_true(all(dim(astats$cos) == c(nrow(mini_nsclc$counts), 3)))
  expect_true(all(dim(astats$llr) == c(nrow(mini_nsclc$counts), 3)))
})

# test anchor selection from stats:
anchors <- choose_anchors_from_stats(counts = mini_nsclc$counts,
                                     neg = Matrix::rowMeans(mini_nsclc$neg),
                                     bg = NULL,
                                     anchorstats = astats, 
                                     cos = NULL, 
                                     llr = NULL, 
                                     n_cells = 500, 
                                     min_cosine = 0.3, 
                                     min_scaled_llr = 0.01, 
                                     insufficient_anchors_thresh = 2) 

testthat::test_that("choose_anchors_from_stats produces correct outputs", {
  expect_true(all(is.element(anchors, c(NA, colnames(astats[[2]])))))
  expect_true(max(table(anchors)) <= 500)
  expect_equal(names(anchors), rownames(mini_nsclc$counts))
})

# test global anchor selection:
anchors <- find_anchor_cells(counts = mini_nsclc$counts,
                             neg = Matrix::rowMeans(mini_nsclc$neg),
                             bg = NULL, 
                             align_genes = TRUE,
                             profiles = ioprofiles[, 1:3], 
                             size = 10, 
                             n_cells = 500, 
                             min_cosine = 0.3, 
                             min_scaled_llr = 0.01, 
                             insufficient_anchors_thresh = 2) 

testthat::test_that("find_anchor_cells produces correct outputs", {
  expect_true(all(is.element(anchors, c(NA, colnames(astats[[2]])))))
  expect_true(max(table(anchors)) <= 500)
  expect_equal(names(anchors), rownames(mini_nsclc$counts))
})

# test anchor selection runs without error and returns NULL if no cells meet criteria:
anchors <- suppressWarnings(find_anchor_cells(counts = mini_nsclc$counts,
                             neg = Matrix::rowMeans(mini_nsclc$neg),
                             bg = NULL, 
                             align_genes = TRUE,
                             profiles = ioprofiles[, 1:3], 
                             size = 10, 
                             n_cells = 500, 
                             min_cosine = 0.8, 
                             min_scaled_llr = 0.01, 
                             insufficient_anchors_thresh = 2))
testthat::test_that("find_anchor_cells produces correct outputs when none selected", {
  expect_null(anchors)
})


# test merge cells:
merge1 <- mergeCells(
  merges = NULL, to_delete = NULL, logliks = sup$logliks)
merge2 <- mergeCells(
  merges =  c("macrophage" = "myeloid", "mDC" = "myeloid","B-cell" = "lymphoid"), 
  to_delete = "endothelial", 
  logliks = sup$logliks)
testthat::test_that("mergecells works when no directions are passed to it", {
  expect_equal(sup$logliks, merge1$logliks, tolerance = 1e-2)
  expect_equal(sup$clust, merge1$clust)
})
testthat::test_that("mergecells works when merges and deletions are asked for", {
  expect_true(all(is.element(colnames(merge2$logliks), c("lymphoid", "myeloid", "fibroblast", "mast"))))
  expect_true(all(is.na(merge2$logliks[sup$logliks[, "endothelial"] == 1, 1])))
  expect_equal(names(merge2$clust), names(sup$clust))
})
