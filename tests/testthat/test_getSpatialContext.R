
# load data ("raw" and "cellannot"):
data("mini_nsclc")


test_that("getNeighborhood expression works under diverse settings", {

  n1 <- getSpatialContext(counts = mini_nsclc$counts, xy = cbind(mini_nsclc$x, mini_nsclc$y), N = 50)
  expect_equal(dim(n1), dim(mini_nsclc$counts))
  
  n2 <- getSpatialContext(counts = mini_nsclc$counts, xy = cbind(mini_nsclc$x, mini_nsclc$y), rad = 0.1, dim_reduce_to = 20)
  expect_equal(dim(n2), c(nrow(mini_nsclc$counts), 20))
})

