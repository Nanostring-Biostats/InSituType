if (FALSE) {
  library(testthat)
  library(MLEcell)
  library(Matrix)
}

data("ioprofiles")
data("mini_nsclc")
mat <- mini_nsclc$counts
bg <- Matrix::rowMeans(mini_nsclc$neg)
genes <- intersect(dimnames(mini_nsclc$counts)[[2]], dimnames(ioprofiles)[[1]])
x <- ioprofiles[genes, 1]
testthat::test_that("negative binomial disrtibution is same as stats package", {
  bgsub <- pmax( sweep( mat , 1 , bg , "-" ) , 0 )
  s <- Matrix::rowSums(bgsub) / sum(x)
  s[s <= 0] <- Matrix::rowSums(mat[s <= 0, , drop = FALSE]) / sum(x)
  result <- MLEcell::lls(as(mat, "dgCMatrix"), s, x, bg, 10)
  names(result) <- rownames(mat)
  yhat <- sweep( s %*% t( x ) , 1 , bg , "+" )
  lls <- stats::dnbinom(x = as.matrix(mat), size = 10, mu = yhat, log = TRUE)
  result_ref <- rowSums(lls)
  expect_true(all.equal(result, result_ref, tolerance = 0.02))
})
