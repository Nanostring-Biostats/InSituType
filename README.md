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

### FAQs and tips:
[https://github.com/Nanostring-Biostats/InSituType/FAQs.md](https://github.com/Nanostring-Biostats/InSituType/blob/main/FAQs.md)

### Reproduction instructions
The full results of the Insitutype manuscript can be reproduced with the code in this repo: https://github.com/Nanostring-Biostats/InSituType-manuscript-analyses

### Installation
```
# Make sure Matrix and irlba are both up to date (otherwise versioning issues cause prcomp_irlba to error out):
# (This is required as of Feb 2024; with any luck these packages will fix their versioning issues soon and this will not be necessary.)
install.packages("Matrix", type = "source")
install.packages("irlba", type = "source")

# Install Insitutype:
devtools::install_github("https://github.com/Nanostring-Biostats/InSituType")
```
Installation should take < 2 mins on a normal desktop computer. 


### Function dependencies:
![image](https://user-images.githubusercontent.com/4357938/200046292-ba3e3453-b201-4776-b5f5-6bf3dfce6ec6.png)
