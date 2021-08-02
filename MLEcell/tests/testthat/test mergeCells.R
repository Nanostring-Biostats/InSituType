

# example logliks:
logliks <- matrix(c(-3,-3,-2,-1,-1,-2), 
                  nrow = 2,
                  dimnames = list(paste0("cell", 1:2), paste0("old_", letters[1:3])))
# define merges:
merges <- c("old_a" = "new1", "old_b" = "new1")

# run:
res <- mergeCells(merges = merges, logliks = logliks)

# confirm it works:
testthat::test_that("new cluster names are right", {
  expect_equal(colnames(res$probs), c("new1", "old_c"))
})

testthat::test_that("new cluster assignments are right", {
  expect_equal(res$clust, c("cell1" = "old_c", "cell2" = "new1"))
})

testthat::test_that("probabilities are right", {
  expect_equal(res$probs[,1], c("cell1" = 0.27, "cell2" = 0.73), tolerance = 2)
})
