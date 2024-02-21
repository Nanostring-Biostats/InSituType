# InSituType 2.0.0

* Enable use in protein datasets via the assay_type argument. This required a major overhaul under the hood, but has little impact on existing RNA workflows. 
* More advanced methods for updating reference profiles via anchor cells, implemented in `updateReferenceProfiles`.
* New function `spatialUpdate` for using alternative data types (e.g. space or immunofluorescence) and the Insitutype likelihood framework to update cell typing results from any method. 
* New functions `getRNAprofiles` and `getProteinParameters`, which serve as user-facing tools for getting profile matrices. 


# InSituType 1.0.0

* License updated
* lldist parallelized with OpenMP

# InSituType 0.99.4

* Re-submission to Bioconductor

# InSituType 0.99.3

* Merge subclustering fix

# InSituType 0.99.2

* Optionally use SingleCellExperiment class

# InSituType 0.99.1

* Added reference to CosMx paper and dataset

# InSituType 0.99.0

* Submission to Bioconductor 3.16

# InSituType 1.1.1

* Updated `flightpath_layout.R` to save the plot in a temp folder in the current work directory

# InSituType 1.1.0

* Fix several places counts matrix was being converted to dense to calculate a statistic
  * Revert conversion to `sparse matrix` of `dense` `mu` matrix and result from `dnbinom`

# InSituType 1.0.0

* Integrated rcpp support for the package
* Added `dnbinom` for `sparse matrices`
* Updated the `unit tests`
* Removed `.o` files in `src` folder

# InSituType 0.1.2

* Updated package dependencies

# InSituType 0.1.1

* Added `lsa`, `SpatialDecon`, `irlba`, `mclust`, `rmarkdown` to the `DESCRIPTION` file
* Fixed a Roxygen example for the R function `geoSketch` where it was trying to use `Ptolemy` and `Giotto` packages that are not being used within the package

# InSituType 0.1.0

* Added a `NEWS.md` file to track changes to the package.
* Added BioConductor package dependencies (notably SpatialDecon and lsa)
* Renamed vignettes to allow for compilation
* Deleted old vignettes (labelled OLD)
