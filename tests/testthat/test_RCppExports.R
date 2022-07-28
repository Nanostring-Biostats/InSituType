if (FALSE) {
  library(testthat)
  library(MLEcell)
}

i <- c(1,3:8)
j <- c(2,9,6:10)
x <- 7 * (1:7)
A <- Matrix::sparseMatrix(i, j, x = x)
x <- 1:7
mu <- Matrix::sparseMatrix(i, j, x = x)
result <- MLEcell::dnbinom_sparse(x=A, mu=mu, size_dnb = 10)
result_n <- stats::dnbinom(x=as.matrix(A), mu=as.matrix(mu), size = 10, log=T)
testthat::test_that("negative binomial disrtibution is same as stats package", {
  expect_true(identical(result_n, as.matrix(result)))
})
