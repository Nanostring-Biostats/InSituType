# load data ("raw" and "cellannot"):
data("ioprofiles")
data("mini_nsclc")

# run unsupervised clustering with several random starts:
unsup <- insitutype(counts = mini_nsclc$counts,
                    neg = rowMeans(mini_nsclc$neg),
                    bg = NULL,
                    init_clust = NULL, n_clusts = 12,
                    fixed_profiles = NULL,
                    nb_size = 10,
                    method = "EM",
                    n_starts = 2,
                    n_benchmark_cells = 500,
                    n_phase1 = 100,
                    n_phase2 = 200,
                    n_phase3 = 300)   
# plot clusters:
scols = c(cellcols, brewer.pal(12, "Set3")[-2], brewer.pal(8, "Set2"))[1:length(unique(unsup$clust))]
names(scols) = unique(unsup$clust)
plot(mini_nsclc$x, mini_nsclc$y, pch = 16, col = scols[unsup$clust])

# view immune oncology cell profiles (in ptolemy package data):
head(ioprofiles)
semi <- insitutype(counts = mini_nsclc$counts,
                   neg = rowMeans(mini_nsclc$neg),
                   bg = NULL,
                   init_clust = NULL, n_clusts = 2,
                   fixed_profiles = ioprofiles,
                   nb_size = 10,
                   method = "EM",
                   n_starts = 2,
                   n_benchmark_cells = 500,
                   n_phase1 = 100,
                   n_phase2 = 200,
                   n_phase3 = 300)   
