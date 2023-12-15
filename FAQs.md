### FAQs and advanced methods

#### Topics
[Choosing the n_clust argument](#choosing-nclust)
[Updating reference profiles](#updating-reference-profiles)



### Workflow overview
The broad Insitutype workflow is as follows:
![image](https://github.com/Nanostring-Biostats/InSituType/assets/4357938/bda58044-0904-4b90-8ffe-b40abdb8e222)


### Choosing nclust
We recommend choosing a slightly generous value of nclust, then using refineClusters to condense the resulting clusters. For example, if you're running semi-supervised cell typing and you expect to find 5 new clusters, set nclust = 8. Or for unsupervised clustering with an expectation of 12 cell types, set nclust = 16. 
It's generally easy to tell when two clusters come from the same cell type: they'll be adjacent in UMAP space, and the flightpath plot will show them frequently confused with each other. 

Final note: Insitutype splits big clusters with higher counts more aggressively than other clusters. For example, in a tumor study, it will subcluster tumor cells many times before it subclusters e.g. fibroblasts. The simplest solution is to increase nclust as needed, then condense the over-clustered cell type as desired. 

### Anchor cell considerations



### Updating reference profiles

Cell typing's biggest challenge is using a reference dataset from a different platform. Platform effects between single cell and spatial platforms can be profound. 
Insitutype has 3 treatments for reference profiles:
1. Use the reference profile matrix as-is
2. Choose anchor cells, then rescale genes based on estimated platform effects. (Less aggressive, only fits gene-level effects.)
3. Choose anchor cells, then refit the reference profiles entirely. (Most aggressive, fits a new value for every gene x cell type.)

We suggest using the below flowchart to choose from among these options:

![image](https://github.com/Nanostring-Biostats/InSituType/assets/4357938/e58f8196-f226-4641-8dfa-bafd9b3dbfae)


### Which genes to use

Insitutype was designed using 1000-plex CosMx data, where we found it most powerful to use all genes in the panel. 
In our new 6000-plex data, it's worth considering using Insitutype on a well-chosen subset of genes. As a rule of thumb, genes should be retained if either of the following applies: 
1. They have solidly above-background expression in the CosMx data
2. They have moderate-to-high expression in at least one reference profile

For typical 6000plex experiments, we speculate that cell typing using somewhere between 3000-5000 genes would be optimal. 


### Interpreting clustering results

Once Insitutype has run, take time to scrutinize the results. You'll need to:
1. Confirm cell types from the reference profiles are correct
2. Interpret new clusters

First, we recommend the following QC plots:

![image](https://github.com/Nanostring-Biostats/InSituType/assets/4357938/aa2c47ba-8c4e-412d-b790-5205ae9739fc)
![image](https://github.com/Nanostring-Biostats/InSituType/assets/4357938/f1f1694c-c0df-41fe-a823-ca34a16d553b)

Example code for generating the above profiles heatmap:
```
pdf("<writehere.pdf>", height = 20, width = 6)
mat <- res$profiles  # ("res" is the insitutype output)
mat <- sweep(mat, 1, pmax(apply(mat, 1 ,max), 0.1), "/")
pheatmap(mat, col = colorRampPalette(c("white", "darkblue"))(100),
         fontsize_row = 5)
dev.off()
```

We have found the below workflows to be effective and efficent:
![image](https://github.com/Nanostring-Biostats/InSituType/assets/4357938/2f7774db-7437-447b-a0e1-f32acb3eafcb)

![image](https://github.com/Nanostring-Biostats/InSituType/assets/4357938/8cc7f603-d8e9-4a83-9479-ebb7d528edd8)



### Targeted subclustering

This is an advanced method. Sometimes it can be hard to subcluster a cell type if manyu of its genes are impacted by contamination from segmentation errors. Immune cells in the context of tumors are a good example.
To subcluster say T-cells in a tumor, you might initially call a single T-cell cluster. Then, considering just these cells and just the genes unlikely to be contaminated in T-cells (genes with high T-cell expression or with low expression in surrounding cell types), run unsupervised Insitutype. 




