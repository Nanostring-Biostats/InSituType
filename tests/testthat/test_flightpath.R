
# load data ("raw" and "cellannot"):
data("ioprofiles")
data("iocolors")
data("mini_nsclc")


# run unsupervised clustering with several random starts:
res <- insitutype(x = mini_nsclc$counts,
                  neg = Matrix::rowMeans(mini_nsclc$neg),
                  bg = NULL,
                  init_clust = NULL, n_clusts = 6,
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
                  max_iters = 2,
                  assay_type="RNA")


# test flightpath_layout
fp <- flightpath_layout(probs = NULL, logliks = res$logliks, profiles = res$profiles)

test_that("flightpath_layout returns correct format", {
  expect_true(all(dim(fp$clustpos) == c(6, 2)))
  expect_true(all(dim(fp$cellpos) == c(nrow(res$logliks), 2)))
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

