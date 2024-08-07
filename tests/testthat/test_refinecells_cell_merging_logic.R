# example logliks:
logliks <- matrix(c(-3, -3, -2, -1, -1, -2), 
                  nrow = 2,
                  dimnames = list(paste0("cell", 1:2), paste0("old_", letters[1:3])))

# define merges:
merges <- c("old_a" = "new1", "old_b" = "new1", "old_c" = "old_c")

# run:
res <- refineClusters(merges = merges, logliks = logliks)

# confirm it works:
test_that("new cluster names are right", {
  expect_equal(colnames(res$logliks), c("new1", "old_c.new"))
})

test_that("new cluster assignments are right", {
  expect_equal(res$clust, c("cell1" = "old_c.new", "cell2" = "new1"))
})

test_that("probabilities are right", {
  expect_equal(res$logliks[, 1], c("cell1" = -2, "cell2" = -1), tolerance = 2)
})
