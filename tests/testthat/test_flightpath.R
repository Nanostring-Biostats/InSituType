
# simulate probs and profiles:
celltypes <- letters[1:6]
probs <- matrix(runif(600), 100)
probs <- sweep(probs, 1, rowSums(probs), "/")
profiles <- matrix(rgamma(300, 1), 50)
colnames(probs) <- colnames(profiles) <- celltypes
res <- list(probs = probs, profiles = profiles, clust = colnames(probs)[apply(probs,1,which.max)])


# test flightpath_layout
fp <- flightpath_layout(res$probs, res$profiles)
test_that("flightpath_layout returns correct format", {
  expect_true(all(dim(fp$clustpos) == c(6,2)))
  expect_true(all(dim(fp$cellpos) == c(100,2)))
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


