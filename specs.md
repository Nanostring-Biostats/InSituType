

#### Specs for insitutypeML:
- Returns a vector of cell type assignments  -- test: test_insitutype.R#L54
- Returns a vector of posterior probabilities / confidence scores   -- test: test_insitutype.R#L55
- Returns a matrix of cell * cell type log-likelihoods  -- test: test_insitutype.R#L56
- Returns a matrix of cell type profiles  -- test: test_insitutype.R#L57


#### Specs for insitutype:
- If run with fixed_profiles and 0 new clusters, produces valid outputs:
  - Returns a vector of cell type assignments  -- test: test_insitutype.R#L82
- Returns a vector of posterior probabilities / confidence scores   -- test: test_insitutype.R#L83
- Returns a matrix of cell * cell type log-likelihoods  -- test: test_insitutype.R#L84
- Returns a matrix of cell type profiles  -- test: test_insitutype.R#L85
- If run with no fixed_profiles (fully unsupervises), produces valid outputs:
  - Returns a vector of cell type assignments  -- test: test_insitutype.R#L112
- Returns a vector of posterior probabilities / confidence scores   -- test: test_insitutype.R#L113
- Returns a matrix of cell * cell type log-likelihoods  -- test: test_insitutype.R#L114
- Returns a matrix of cell type profiles  -- test: test_insitutype.R#L115
- If unsupervised clustering is run with initial clusters specified, produces valid outputs:
  - Returns a vector of cell type assignments  -- test: test_insitutype.R#L144
- Returns a vector of posterior probabilities / confidence scores   -- test: test_insitutype.R#L145
- Returns a matrix of cell * cell type log-likelihoods  -- test: test_insitutype.R#L146
- Returns a matrix of cell type profiles  -- test: test_insitutype.R#L147
- The clusters returned have the same names as the initial clusters  -- test: test_insitutype.R#L148
- If semi-supervised clustering is run with initial clusters specified, produces valid outputs:
  - Returns a vector of cell type assignments  -- test: test_insitutype.R#L177
- Returns a vector of posterior probabilities / confidence scores   -- test: test_insitutype.R#L178
- Returns a matrix of cell * cell type log-likelihoods  -- test: test_insitutype.R#L179
- Returns a matrix of cell type profiles  -- test: test_insitutype.R#L180


#### Specs for updateReferenceProfiles:
- Returns a matrix of new profiles  -- test: test_insitutype.R#L316
- Returns a vector of anchor assignments  -- test: test_insitutype.R#L317


#### Specs for refineClusters:
- Merging operations happen correctly  -- test: test_insitutype.R#L196
- Cell names are preserved  -- test: test_insitutype.R#L196
- Makes no changes if none are requested  -- test: test_insitutype.R#L302
- Merging operations happen correctly if merges and deletions are asked for  -- test: test_insitutype.R#L307
- Merging operations happen correctly if merges are asked for  -- test: test_refinecells_cell_merging_logic.R#L14,18,22


#### Specs for chooseClusterNumber:
- Returns a single value for "best cluster number"  -- test: test_insitutype.R#L219
- Reports the cluster numbers considered  -- test: test_insitutype.R#L220
- Reports the log likelihood from each cluster number #221
- Reports the AIC from each cluster number #222
- Reports the BIC from each cluster number #223


#### Specs for get_anchor_stats
- Returns a matrix of cosine distances  -- test: test_insitutype.R#L236
- Returns a matrix of log likelihood ratios  -- test: test_insitutype.R#L237

#### Specs for choose_anchors_from_stats
- Assigns values consistent with the cell type names of the inputs  -- test: test_insitutype.R#L253
- Assigns no more than the specified number of anchors per cell type  -- test: test_insitutype.R#L253
- The anchors vector aligns to the rows of the counts matrix (cells)  -- test: test_insitutype.R#L254

#### Specs for find_anchor_cells
- Assigns values consistent with the cell type names of the inputs  -- test: test_insitutype.R#L271
- Assigns no more than the specified number of anchors per cell type  -- test: test_insitutype.R#L272
- The anchors vector aligns to the rows of the counts matrix (cells)  -- test: test_insitutype.R#L273
- Returns NULL if no cells meet anchor criteria  -- test: test_insitutype.R#L288


#### Specs for flightpath_layout
- Returns correctly formatted results:
  - Cluster positions are in a 2-column matrix  -- test: test_flightpath.R#L34
- Cell positions are in a 2-column matrix  -- test: test_flightpath.R#L35
- There are no missing cluster positions  -- test: test_flightpath.R#L36
- There are no missing cell positions  -- test: test_flightpath.R#L37

#### Specs for flightpath_plot
- when passed a result from flightpath_layout, flightpath_plot returns a ggplot object  -- test: test_flightpath.R#L43
- when passed an insitutype results, flightpath_plot returns a ggplot object  -- test: test_flightpath.R#L51
- when asked to show meanConfidence, flightpath_plot returns a ggplot object  -- test: test_flightpath.R#L57


#### Specs for fastCohorting
- Returns a vector of cohort assignments  -- test: test_insitutype.R#L325
- Returns the specified number of unique cohorts  -- test: test_insitutype.R#L325
