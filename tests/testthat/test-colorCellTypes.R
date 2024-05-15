# create mock cell type abundances:
data("iocolors")

set.seed(0)
cells <- sample(c(letters[1:10], names(iocolors)[1:4]), 100, replace = TRUE)
tab <- table(cells)

# run using just names: 
cols_names <- colorCellTypes(names = names(tab), freqs = NULL, init_colors = NULL, max_sum_rgb = 600, palette = "brewers") 

# run using abundance info: 
cols_freqs <- colorCellTypes(names = NULL, freqs = tab, init_colors = NULL, max_sum_rgb = 600) 


# test that pre-specified colors are used
cols_init <- colorCellTypes(names = NULL, freqs = tab, init_colors = iocolors, max_sum_rgb = 600) 

test_that("pre-specified colors are used", {
  sharedcells <- intersect(names(tab), names(iocolors))
  expect_true(all.equal(cols_init[sharedcells], iocolors[sharedcells]))
})



# test that legal colors are returned in all cases, with names matching the cell names:
test_that("test that results returned by flagLowGenes have the right formats", {
  expect_error(plot(seq_along(tab), col = cols_names), NA) # "NA" means expecting no error
  expect_error(plot(seq_along(tab), col = cols_freqs), NA) # "NA" means expecting no error
  expect_error(plot(seq_along(tab), col = cols_init), NA) # "NA" means expecting no error
  expect_equal(length(intersect(names(cols_names), names(tab))), length(names(tab)))
  expect_equal(length(intersect(names(cols_freqs), names(tab))), length(names(tab)))
  expect_equal(length(intersect(names(cols_init), names(tab))), length(names(tab)))
  
})

# test that it works if prespecified colors have no overlap:
test_that("correct results even if prespecified colors have no overlap", {
  cols_bad_init <- colorCellTypes(names = NULL, freqs = tab, init_colors = c(no = "red", nope = "blue"), max_sum_rgb = 600) 
  expect_error(plot(seq_along(tab), col = cols_bad_init), NA) # "NA" means expecting no error
  expect_equal(length(intersect(names(cols_bad_init), names(tab))), length(names(tab)))
})


# test that all 3 paletted work:
test_that("all 3 paletted work", {
  cols_tab20 <- colorCellTypes(names = NULL, freqs = tab, init_colors = NULL, max_sum_rgb = 600, palette = "tableau20")
  cols_brew <- colorCellTypes(names = NULL, freqs = tab, init_colors = NULL, max_sum_rgb = 600, palette = "brewers") 
  cols_earth <- colorCellTypes(names = NULL, freqs = tab, init_colors = NULL, max_sum_rgb = 600, palette = "earthplus") 
  expect_error(plot(seq_along(tab), col = cols_tab20), NA) # "NA" means expecting no error
  expect_error(plot(seq_along(tab), col = cols_brew), NA) # "NA" means expecting no error
  expect_error(plot(seq_along(tab), col = cols_earth), NA) # "NA" means expecting no error
})

