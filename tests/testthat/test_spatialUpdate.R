data("ioprofiles")
data("iocolors")
data("mini_nsclc")


initclust <- sample(c("a","b","c"), nrow(mini_nsclc$counts), replace = TRUE)

updatedclust <- spatialUpdate(celltype = initclust, 
                              counts = mini_nsclc$counts,
                              neg = Matrix::rowMeans(mini_nsclc$neg),
                              cohort = NULL, altdata = NULL, 
                              xy = cbind(mini_nsclc$x, mini_nsclc$y), 
                              tissue = NULL,
                              nb_size = 10, assay_type = "rna")
test_that("spatialUpdate worked", {
  expect_true(all(is.element(c( "clust","prob","profiles","sds","logliks","logliks_from_lost_celltypes"), names(updatedclust))))
})


