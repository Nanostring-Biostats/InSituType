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


### Targeted subclustering


### Choosing nclust
### Choosing nclust


