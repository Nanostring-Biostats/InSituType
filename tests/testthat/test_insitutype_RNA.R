# load data ("raw" and "cellannot"):
data("ioprofiles")
data("iocolors")
data("mini_nsclc")

set.seed(0)


# test V1.0 syntax still works for back-compatibility:
test_that("old syntax semi-sup still runs", {
  sharedgenes <- intersect(colnames(mini_nsclc$counts), rownames(ioprofiles))
  semisup <- insitutype(
    x = mini_nsclc$counts[, sharedgenes],
    neg = Matrix::rowMeans(mini_nsclc$neg),
    n_clusts = c(5, 6),
    reference_profiles = ioprofiles[sharedgenes, 1:3],
    update_reference_profiles = FALSE,
    max_iters = 3
  ) 
  
  expect_true(all(is.element(c("clust","prob","profiles","logliks","logliks_from_lost_celltypes"), names(semisup))))
})


# test nbclust semi-sup
sharedgenes <- intersect(colnames(mini_nsclc$counts), rownames(ioprofiles))
nbres <- nbclust(counts = mini_nsclc$counts[, sharedgenes],
                 neg =  Matrix::rowMeans(mini_nsclc$neg), 
                 bg = NULL,
                 assay_type = "rna",
                 fixed_profiles = ioprofiles[sharedgenes, 1:3],
                 init_profiles = NULL, 
                 init_clust = rep(c("a", "b"), nrow(mini_nsclc$counts) / 2),
                 nb_size = 10,
                 cohort = rep("a", nrow(mini_nsclc$counts)),
                 pct_drop = 1/10000,
                 min_prob_increase = 0.05, 
                 max_iters = 3, 
                 logresults = FALSE)
test_that("semi-sup nbclust preserves fixedprofiles", {
  expect_true(all(abs(diag(cor(nbres$profiles[, colnames(ioprofiles)[1:3]], ioprofiles[sharedgenes, ]))) == 1))
})


# test supervised cell typing using direct loglik calcs:
sup <- insitutypeML(x = mini_nsclc$counts,
                    neg = Matrix::rowMeans(mini_nsclc$neg),
                    bg = NULL,
                    assay_type = "rna",
                    cohort = rep(c("a", "b"), each = nrow(mini_nsclc$counts) / 2),
                    reference_profiles = ioprofiles[, 1:6],
                    nb_size = 10,
                    align_genes = TRUE)
  
test_that("supervised cell typing produces correct outputs", {
  expect_true(all(is.element(c("clust", "prob", "logliks", "profiles"), names(sup))))
  expect_true(is.vector(sup$clust))
  expect_true(is.vector(sup$prob))
  expect_true(is.matrix(sup$logliks))
  expect_true(is.matrix(sup$profiles))
})

# test semi-supervised with 0 new clusts:
semi <- insitutype(x = mini_nsclc$counts,
                  neg = Matrix::rowMeans(mini_nsclc$neg),
                  bg = NULL,
                  assay_type = "rna",
                  init_clust = NULL,
                  n_clusts = 0,
                  anchors = NULL,
                  reference_profiles = ioprofiles[, 1:6],
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
                  n_anchor_cells = 20,
                  min_anchor_cosine = 0.3,
                  min_anchor_llr = 0.01,
                  insufficient_anchors_thresh = 2,
                  refinement = FALSE, 
                  rescale = FALSE, 
                  refit = TRUE)

test_that("semiservised cell typing with n_clusts = 0 produces correct outputs", {
  expect_true(all(is.element(c("clust", "prob", "logliks", "profiles"), names(semi))))
  expect_true(is.vector(semi$clust))
  expect_true(is.vector(semi$prob))
  expect_true(is.matrix(semi$logliks))
  expect_true(is.matrix(semi$profiles))
})


# run unsupervised clustering with several random starts:
unsup <- insitutype(x = mini_nsclc$counts,
                    neg = Matrix::rowMeans(mini_nsclc$neg),
                    bg = NULL,
                    assay_type = "rna",
                    init_clust = NULL, 
                    n_clusts = 2:5,
                    reference_profiles = NULL,
                    reference_sds = NULL,
                    anchors = NULL,
                    cohort = rep(c("a", "b"), each = nrow(mini_nsclc$counts) / 2),
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
                    max_iters = 2,
                    refinement = FALSE, 
                    rescale = FALSE, 
                    refit = TRUE)   

test_that("unsupervised cell typing produces correct outputs", {
  expect_true(all(is.element(c("clust", "prob", "logliks", "profiles"), names(unsup))))
  expect_true(is.vector(unsup$clust))
  expect_true(is.vector(unsup$prob))
  expect_true(is.matrix(unsup$logliks))
  expect_true(is.matrix(unsup$profiles))
})




# run unsupervised clustering with init_clust specified:
init_clust <- rep(c("name1", "xxx", "ooo"), each = nrow(mini_nsclc$counts) / 3)[seq_len(nrow(mini_nsclc$counts))]
unsup <- insitutype(x = mini_nsclc$counts,
                    neg = Matrix::rowMeans(mini_nsclc$neg),
                    bg = 0.03,
                    assay_type = "rna",
                    init_clust = init_clust, 
                    n_clusts = 6,
                    reference_profiles = NULL,
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
                    max_iters = 4,
                    refinement = FALSE, 
                    rescale = FALSE, 
                    refit = TRUE)   


test_that("unsupervised cell typing using init_clust produces correct outputs", {
  expect_true(all(is.element(c("clust", "prob", "logliks", "profiles"), names(unsup))))
  expect_true(is.vector(unsup$clust))
  expect_true(is.vector(unsup$prob))
  expect_true(is.matrix(unsup$logliks))
  expect_true(is.matrix(unsup$profiles))
  expect_true(all(sort(unique(unsup$clust)) == sort(unique(init_clust))))
})


# semi-supervised using the immune oncology cell profiles (in ptolemy package data):
semi <- insitutype(x = mini_nsclc$counts,
                   neg = Matrix::rowMeans(mini_nsclc$neg),
                   bg = NULL,
                   anchors = NULL,
                   assay_type = "rna",
                   init_clust = NULL, n_clusts = 2,
                   reference_profiles = ioprofiles[, 1:3],
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
                   sketchingdata = NULL, 
                   insufficient_anchors_thresh = 2,
                   refinement = FALSE, 
                   rescale = FALSE, 
                   refit = TRUE)   


test_that("semi-supervised cell typing using init_clust produces correct outputs", {
  expect_true(all(is.element(c("clust", "prob", "logliks", "profiles"), names(semi))))
  expect_true(is.vector(semi$clust))
  expect_true(is.vector(semi$prob))
  expect_true(is.matrix(semi$logliks))
  expect_true(is.matrix(semi$profiles))
})


# test merge cells with multi-sample clustering:
merge2 <- refineClusters(
  merges =  c("B-cell" = "lymphoid", "a" = "cancer", "b" = "cancer"), 
  to_delete = c("endothelial"), 
  subcluster = list("fibroblast" = 2),
  logliks = semi$logliks, 
  counts = mini_nsclc$counts,
  neg = Matrix::rowMeans(mini_nsclc$neg),
  bg = NULL, 
  cohort = NULL)
test_that("refineClusters works when merges and deletions are asked for", {
  expect_true(all(is.element(colnames(merge2$logliks), c("lymphoid", "cancer", "fibroblast_1", "fibroblast_2"))))
  expect_equal(names(merge2$clust), names(semi$clust))
})



# test chooseclusternumber:
res <- chooseClusterNumber(counts = mini_nsclc$counts,
                           neg = Matrix::rowMeans(mini_nsclc$neg),
                           bg = NULL,
                           anchors = NULL,
                           assay_type = "rna",
                           init_clust = NULL, n_clusts = 2:3,
                           fixed_profiles = ioprofiles[, 1:3],
                           max_iters = 4, 
                           subset_size = 1000, 
                           align_genes = TRUE, 
                           plotresults = FALSE, 
                           nb_size = 10, 
                           pct_drop = 0.005, 
                           min_prob_increase = 0.05) 

test_that("chooseClusterNumber produces correct outputs", {
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
                           sds=NULL,
                           assay_type = "rna",
                           align_genes = TRUE,
                           profiles = ioprofiles[, 1:3], 
                           size = 10, 
                           min_cosine = 0.3) 

test_that("get_anchor_stats produces correct outputs", {
  expect_true(all(dim(astats$cos) == c(nrow(mini_nsclc$counts), 3)))
  expect_true(all(dim(astats$llr) == c(nrow(mini_nsclc$counts), 3)))
})

# test anchor selection from stats:
anchors <- choose_anchors_from_stats(counts = mini_nsclc$counts,
                                     neg = Matrix::rowMeans(mini_nsclc$neg),
                                     bg = NULL,
                                     assay_type = "rna",
                                     anchorstats = astats, 
                                     cos = NULL, 
                                     llr = NULL, 
                                     n_cells = 500, 
                                     min_cosine = 0.3, 
                                     min_scaled_llr = 0.01, 
                                     insufficient_anchors_thresh = 2) 

test_that("choose_anchors_from_stats produces correct outputs", {
  expect_true(all(is.element(anchors, c(NA, colnames(astats[[2]])))))
  expect_true(max(table(anchors)) <= 500)
  expect_equal(names(anchors), rownames(mini_nsclc$counts))
})

# test global anchor selection:
anchors <- find_anchor_cells(counts = mini_nsclc$counts,
                             neg = Matrix::rowMeans(mini_nsclc$neg),
                             bg = NULL, 
                             align_genes = TRUE,
                             assay_type = "rna",
                             profiles = ioprofiles[, 1:3], 
                             sds=NULL,
                             size = 10, 
                             n_cells = 500, 
                             min_cosine = 0.3, 
                             min_scaled_llr = 0.01, 
                             insufficient_anchors_thresh = 2) 

test_that("find_anchor_cells produces correct outputs", {
  expect_true(all(is.element(anchors, c(NA, colnames(astats[[2]])))))
  expect_true(max(table(anchors)) <= 500)
  expect_equal(names(anchors), rownames(mini_nsclc$counts))
})

# test anchor selection runs without error and returns NULL if no cells meet criteria:
anchors <- suppressWarnings(find_anchor_cells(counts = mini_nsclc$counts,
                             neg = Matrix::rowMeans(mini_nsclc$neg),
                             bg = NULL, 
                             align_genes = TRUE,
                             assay_type = "rna",
                             profiles = ioprofiles[, 1:3], 
                             sds=NULL,
                             size = 10, 
                             n_cells = 500, 
                             min_cosine = 0.8, 
                             min_scaled_llr = 0.01, 
                             insufficient_anchors_thresh = 2))
test_that("find_anchor_cells produces correct outputs when none selected", {
  expect_null(anchors)
})


# test refineClusters:
merge1 <- refineClusters(
  merges = NULL, to_delete = NULL, subcluster = NULL, logliks = sup$logliks)
test_that("refineClusters works when no directions are passed to it", {
  expect_equal(sup$logliks, merge1$logliks, tolerance = 1e-2)
  expect_equal(sup$clust, merge1$clust)
})

merge2 <- refineClusters(
  merges =  c("macrophage" = "myeloid", "mDC" = "myeloid", "B-cell" = "lymphoid"), 
  to_delete = "endothelial", 
  subcluster = list("fibroblast" = 2),
  logliks = sup$logliks,
  counts = mini_nsclc$counts,
  neg = Matrix::rowMeans(mini_nsclc$neg))
test_that("refineClusters works when merges and deletions are asked for", {
  expect_true(all(is.element(colnames(merge2$logliks), c("lymphoid", "myeloid", "fibroblast_1", "fibroblast_2", "mast"))))
})


# test updateReferenceProfiles:
rescaled <- updateReferenceProfiles(reference_profiles = ioprofiles, 
                                    reference_sds = NULL,
                                    assay_type = "rna",
                                    counts = mini_nsclc$counts,
                                    neg = Matrix::rowMeans(mini_nsclc$neg))

test_that("rescaleProfiles works as intended", {
  expect_true(is.matrix(rescaled$updated_profiles))
  expect_true(is.vector(rescaled$anchors))
})


# test fastCohorting:
simdat <- matrix(rnorm(10000), nrow = 1000)
cres <- fastCohorting(simdat, n_cohorts = 10)
test_that("fastCohorting works as intended", {
  expect_true(is.vector(cres))
  expect_true(length(unique(cres)) == 10)
})
