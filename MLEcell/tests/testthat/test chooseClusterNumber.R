library(testthat)
library(Ptolemy)
data("mini_tma")

counts <- t(as.matrix(mini_tma@expression$rna$raw[1:100, ]))
neg <- Matrix::colMeans(mini_tma@expression$neg$raw)

test_that("chooseClusterNumber works with no fixed profiles", {
  expect_error(
    res <- chooseClusterNumber(counts = counts, neg = neg,
                               n_clusts = 2:4,
                               n_iters = 4, subset_size = 500), NA)

})


test_that("chooseClusterNumber works with fixed profiles", {
  data(ioprofiles)
  expect_error(
    res2 <- chooseClusterNumber(counts = counts, neg = neg,
                                n_clusts = 2:4,
                                n_iters = 4, subset_size = 100,
                                fixed_profiles = ioprofiles[, 1:4],
                                plotresults = TRUE) , NA)
})


# test that a good choice of number of clusters is made:
test_that("chooseClusterNumber makes a sane choice", {
  res3 <- chooseClusterNumber(counts = counts, neg = neg,
                             n_clusts = c(2,6,20),
                             n_iters = 4, subset_size = 500)
  expect_true(res$best_clust_number == 20)
})

