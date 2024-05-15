data("ioprofiles")
data("iocolors")
data("mini_nsclc")


initclust <- sample(c("a","b","c"), nrow(mini_nsclc$counts), replace = TRUE)

test_that("getRNAprofiles worked", {
  temp <- getRNAprofiles(x = mini_nsclc$counts, neg = 0, clust = initclust)
  expect_identical(rownames(temp), colnames(mini_nsclc$counts))
  expect_identical(colnames(temp)[order(colnames(temp))], unique(initclust)[order(unique(initclust))])
})


test_that("getproteinparameters worked", {
  temp <- getProteinParameters(x = mini_nsclc$counts, clust = initclust)
  expect_identical(rownames(temp$profiles), colnames(mini_nsclc$counts))
  expect_identical(rownames(temp$sds), colnames(mini_nsclc$counts))
  expect_identical(colnames(temp$profiles)[order(colnames(temp$profiles))], unique(initclust)[order(unique(initclust))])
  expect_identical(colnames(temp$sds)[order(colnames(temp$sds))], unique(initclust)[order(unique(initclust))])
})


