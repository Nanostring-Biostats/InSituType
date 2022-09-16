



#### Reqs for insitutype:
Insitutype performs unsupervised clustering, or semi-supervised clustering if
provided with reference profiles. It uses an Expectation_maximization (EM) algorithm based on a negbinom 
distribution. Insitutype coordinates calls to nbclust(), which runs the EM algorithm.

##### Inputs:
- an expression matrix (cells * genes)
- a vector of mean negprobe values
- for semi-supervised learning, a matrix of reference profiles
- additional arguments for finer control

##### Outputs:
A list, with the following elements:
\enumerate{
\item clust: a vector given cells' cluster assignments
\item prob: a vector giving the confidence in each cell's cluster
\item logliks: Matrix of cells' log-likelihoods under each cluster. Cells in rows, clusters in columns.
\item profiles: a matrix of cluster-specific expression profiles
\item anchors: from semi-supervised clustering: a vector giving the identifies and cell types of anchor cells
}


#### Reqs for insitutypeML:
Insitutype performs supervised cell typing using a Bayes classifier based on a negbinom distribution. 

##### Inputs:
- an expression matrix (cells * genes)
- a vector of mean negprobe values
- for semi-supervised learning, a matrix of reference profiles
- additional arguments for finer control

##### Outputs:
A list, with the following elements:
\enumerate{
\item clust: a vector given cells' cluster assignments
\item prob: a vector giving the confidence in each cell's cluster
\item logliks: Matrix of cells' log-likelihoods under each cluster. Cells in rows, clusters in columns.
\item profiles: a matrix of cluster-specific expression profiles
}


#### Reqs for updateReferenceProfiles
Update reference profiles from alternative platforms to better fit the spatial platform. 
Uses pre-specified anchor cells, or if no anchors are specified, by first choosing anchor cells.

##### Inputs:
- reference profiles
- spatial data: counts matrix, negmean values
- additional arguments for finer control

##### Outputs:
- An updated reference matrix
- A vector storing the anchor cells used

#### Reqs for refineClusters
A function for refining the output of insitutype and insitutypeML. 
Can delete clusters, merge/rename clusters, or sub-cluster clusters. 

##### Inputs:
- Results from an insitutyle/insitutypeML run
- If subclustering further, counts data

##### Outputs:
A list in the format of insitutype results with updated cluster assignments. 



#### Reqs for chooseClusterNumber
A function to run insituytpe across a range of cluster numbers and identify the best fit

##### Inputs:
- The standard insitutype inputs
- A range of cluster numbers

##### Outputs:
- A suggested cluster number, plus metrics for comparing cluster numbers.




#### Reqs for get_anchor_stats
Function to calculate the summary stats used by anchor cell selection. 
Results are meant to be fed to choose_anchors_from_stats().

##### Inputs:
- The same expression data used by insitutype.
- Reference profiles

##### Outputs:
- A matrix of cosine distances of cells * cell types
- A matrix of log likelihood ratio scores for cells * cell types



#### Reqs for choose_anchors_from_stats
Chooses anchor cells given cosine distances and log likelihood ratio scores 
output by get_anchor_stats. 

##### Inputs:
- A matrix of cosine distances of cells * cell types
- A matrix of log likelihood ratio scores for cells * cell types

##### Outputs:
A vector of anchor assignments. 



#### Reqs for find_anchor_cells
Complete anchor cell selection workflow. Calls get_anchor_stats and choose_anchors_from_stats.

##### Inputs:
- The same expression data used by insitutype.
- Reference profiles

##### Outputs:
A vector of anchor assignments. 



#### Reqs for flightpath_layout
A function to define the layout for a flightpath plot. Uses UMAP to place cluster centroids,
then places cells based on their posterior probabilities of belonging to each centroid.

##### Inputs:
- A matrix of cell * cluster log-likelihoods (output by insitutype)
- A matrix of cluster profiles

##### Outputs:
- xy placements for cluster centroids
- xy placements for individual cells



#### Reqs for flightpath_plot
Makes a ggplot object holding a flightpath plot. Uses UMAP to place cluster centroids,
then places cells based on their posterior probabilities of belonging to each centroid.

##### Inputs:
- Path 1: input an insitutype/insitutypeML result, and it will call flightpath_layout()
- Path 2: input a flightpath_layout result. 

##### Outputs:
A ggplot object


#### Reqs for fastCohorting
Quickly clusters data from alternative sources like immunofluorescence and spatial context. 

##### Inputs:
- A matrix holding alternative data (cells * variables)
- Arguments for finer control

##### Output:
A vector giving each cell's cohort assignment. 
