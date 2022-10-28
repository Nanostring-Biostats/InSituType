# InSituType
 An R package for performing cell typing in SMI and other single cell data

**Manuscript**: https://www.biorxiv.org/content/10.1101/2022.10.19.512902v1.abstract

**Citing Insitutype**: Danaher P, Zhao E, Yang Z, Ross D, Gregory M, Reitz Z, Kim TK, Baxter S, Jackson S, He S, Henderson DA. Insitutype: likelihood-based cell typing for single cell spatial transcriptomics. bioRxiv. 2022 Jan 1.

### System requirements
- R (>= 3.5.0)
- UNIX, Mac or Windows
- Rcpp library (>= 1.0.9)
- see DESCRIPTION for full dependencies

### Demo
See the "vignettes" folder. Vignettes should run in <5 minutes. 

### Instructions for use
Run "insitutype" for unsupervised or semi-supervised clustering. Run "insitutypeML" for supervised cell typing. See the vignettes for example workflows. 

### Reproduction instructions
The full results of the Insitutype manuscript can be reproduced with the code in this repo: https://github.com/Nanostring-Biostats/InSituType-manuscript-analyses

### Installation
```
devtools::install_github("https://github.com/Nanostring-Biostats/InSituType")
```
Installation should take < 2 mins on a normal desktop computer. 


### Organization of package:
![image](https://user-images.githubusercontent.com/4357938/144138602-595a3686-164a-4127-a35d-eb97f4331e4b.png)
