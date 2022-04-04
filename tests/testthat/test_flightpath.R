
# load data ("raw" and "cellannot"):
data("ioprofiles")
data("iocolors")
data("mini_nsclc")


# run unsupervised clustering with several random starts:
res <- insitutype(counts = mini_nsclc$counts,
                    neg = Matrix::rowMeans(mini_nsclc$neg),
                    bg = NULL,
                    tissue =  sample(letters[1:2], nrow(mini_nsclc$counts), replace = TRUE),
                    init_clust = NULL, n_clusts = 2,
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



if (FALSE) {
  # simulate probs and profiles:
  celltypes <- letters[1:6]
  probs <- matrix(runif(600), 100)
  probs <- sweep(probs, 1, rowSums(probs), "/")
  profiles <- matrix(rgamma(300, 1), 50)
  colnames(probs) <- colnames(profiles) <- celltypes
  logliks <- MLEcell:::probs2logliks(probs)
  res <- list(logliks = logliks, profiles = profiles, clust = colnames(logliks)[apply(logliks,1,which.max)])
}

# test flightpath_layout
fp <- flightpath_layout(probs = NULL, logliks = res$logliks, profiles = res$profiles)

#plot(fp$cellpos)
#text(fp$clustpos[,1], fp$clustpos[,2], rownames(fp$clustpos), col = "red")

test_that("flightpath_layout returns correct format", {
  expect_true(all(dim(fp$clustpos) == c(4,2)))
  expect_true(all(dim(fp$cellpos) == c(nrow(res$logliks),2)))
  expect_true(all(!is.na(fp$clustpos)))
  expect_true(all(!is.na(fp$cellpos)))
})


# test flightpath_plot from flightpath results
p <- flightpath_plot(flightpath_result = fp)
test_that("flightpath_plot returns a ggplot object", {
  expect_true(any(grepl("gg", class(p))))
})


# test flightpath_plot from insitutype results
p <- flightpath_plot(insitutype_result = res)
test_that("flightpath_plot returns a ggplot object", {
  expect_true(any(grepl("gg", class(p))))
})

# test flightpath_plot showing meanconfidence
p <- flightpath_plot(insitutype_result = res, showclusterconfidence = TRUE)
test_that("flightpath_plot returns a ggplot object when showclusterconfidence = TRUE", {
  expect_true(any(grepl("gg", class(p))))
})


