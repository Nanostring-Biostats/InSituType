data("ioprofiles")
data("mini_nsclc")
bg <- Matrix::rowMeans(mini_nsclc$neg)
genes <- intersect(dimnames(mini_nsclc$counts)[[2]], dimnames(ioprofiles)[[1]])
mat <- mini_nsclc$counts[, genes]
x <- ioprofiles[genes, 1, drop = FALSE]

test_that("Rcpp calculation is same as stats package for RNA data type", {
  bgsub <- pmax(sweep(mat, 1, bg, "-"), 0)
  s <- Matrix::rowSums(bgsub) / sum(x)
  s[s <= 0] <- Matrix::rowSums(mat[s <= 0, , drop = FALSE]) / sum(x)
  result <- lldist(mat = as(mat, "dgCMatrix"),
                               assay_type = "rna",
                               x = x,
                               bg=bg, 
                               size = 10)
  names(result) <- rownames(mat)
  yhat <- sweep(s %*% t(x), 1, bg, "+")
  lls <- stats::dnbinom(x = as.matrix(mat), size = 10, mu = yhat, log = TRUE)
  result_ref <- round(rowSums(lls), digits=2)
  expect_true(all.equal(result[,1], result_ref))
})


data("tonsil_protein")
data("tonsil_reference_profile")
bg <- Matrix::rowMeans(tonsil_protein$neg)
proteins <- intersect(dimnames(tonsil_protein$counts)[[2]], dimnames(tonsil_reference_profile$mean.ref.profile)[[1]])
mat <- tonsil_protein$counts[, proteins]
x <- tonsil_reference_profile$mean.ref.profile[proteins, 1, drop = FALSE]
xsd <- tonsil_reference_profile$SDs.ref.profile[proteins, 1, drop = FALSE]


test_that("Rcpp calculation is same as stats package for protein data type", {
  bgsub <- pmax(sweep(mat, 1, bg, "-"), 0)
  s <- Matrix::rowSums(bgsub) / sum(x)
  s[s <= 0] <- Matrix::rowSums(mat[s <= 0, , drop = FALSE]) / sum(x)
  result <- lldist(mat = as.matrix(mat),
                               assay_type = "Protein",
                               x = x,
                               xsd = xsd,
                               bg=bg, 
                               size = 10)
  names(result) <- rownames(mat)
  
  yhat <- s %*% t(x)
  ysd <- s %*% t(xsd)

  lls <- stats::dnorm(x = as.matrix(mat), sd = ysd, mean = yhat, log = TRUE)
  
  result_ref <- round(rowSums(lls), digits=2)
  expect_true(all.equal(result[,1], result_ref))
})
