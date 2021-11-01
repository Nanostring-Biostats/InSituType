# load data ("raw" and "cellannot"):
data(ioprofiles)
data("mini_smi")

# run unsupervised clustering with several random starts:
unsup <- insitutype(counts = raw,
                    neg = runif(nrow(raw), 0, 0.1),
                     bg = NULL,
                     init_clust = NULL, n_clusts = 12,
                     fixed_profiles = NULL,
                     nb_size = 10,
                     n_iters = 10,  # this is not enough
                     method = "EM", shrinkage = 0.5,
                     subset_size = 1000,  # smaller than ideal
                     n_starts = 4,
                     n_benchmark_cells = 500,
                     n_final_iters = 5)   # this is not enough
# plot clusters:
scols = c(cellcols, brewer.pal(12, "Set3")[-2], brewer.pal(8, "Set2"))[1:length(unique(unsup$clust))]
names(scols) = unique(unsup$clust)
plot(cellannot$x, cellannot$y, pch = 16, col = scols[unsup$clust])
#'
# view immune oncology cell profiles (in ptolemy package data):
head(ioprofiles)
usegenes = intersect(rownames(ioprofiles), colnames(raw))
semi <- cellEMClust(counts = raw[, usegenes],
                    bg = bg.predicted,
                    init_clust = NULL, n_clusts = 6,
                    fixed_profiles = ioprofiles[usegenes, ],
                    nb_size = 10,
                    n_iters = 10,  # this is not enough
                    method = "EM",
                    shrinkage = 0.5,
                    subset_size = 1000,
                    n_starts = 4,
                    n_benchmark_cells = 500,
                    n_final_iters = 5)  # this is not enough