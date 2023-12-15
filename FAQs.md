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
1. Use as-is
2. Choose anchor cells, then rescale genes based on estimated platform effects. (Less aggressive, only fits gene-level effects.)
3. Choose anchor cells, then refit the reference profiles entirely. (Most aggressive, fits a new value for every gene x cell type.)


![image](https://github.com/Nanostring-Biostats/InSituType/assets/4357938/e58f8196-f226-4641-8dfa-bafd9b3dbfae)


### Targeted subclustering


### Which genes to use

Insitutype was designed using 1000-plex CosMx data, where we found it most powerful to use all genes in the panel. 
In our new 6000-plex data, it's worth considering using Insitutype on a well-chosen subset of genes. As a rule of thumb, genes should be retained if either of the following applies: 
1. They have solidly above-background expression in the CosMx data
2. They have moderate-to-high expression in at least one reference profile

For typical 6000plex experiments, we speculate that cell typing using somewhere between 3000-5000 genes would be optimal. 


### Choosing nclust


