# General
This repository contains the code used for my thesis for MSc in Statistics and Data Science. The title of the thesis is the following: Penalized reduced rank regression for multidimensional survival data: new estimation approaches and simulation study.

# Folders information
There are 3 folders:
## SimultionStudy1
This folder contains the R markdown file used for Simulation Study 1. In the workspace, the results are saved, except the datasets, since those can be easily reproduced fast by running the specific chunks of the Markdown file. The plots generated for simulation study 1 are also in this folder.
## SimultionStudy2
This folder contains the R markdown file used for Simulation Study 2. In the workspace, the results are saved, except the datasets, since those can be easily reproduced fast by running the specific chunks of the Markdown file. The relevant.Rds files are also included in the folder as well as the plots that were generated for the thesis.
## LASSO-penalized_survRRR_ADMM_implementation
In this folder, the R markdown file which was used for the comparisons of ADMM-based survRRR, glmnet-based survRRR, and unpenalized survRRR algorithms is presented. The ADMM algorithm is presented in a chick of the markdown file while the glmnet-based algorithm is in pen_survRRR.R file and was taken from the study of Sluiskes et al. (2024). 

# References:
1. Sluiskes, M. H., Putter, H., Beekman, M., Goeman, J. J., & Rodríguez-Girondo, M. (2024). Penalized reduced rank regression for multi-outcome survival data supports a common metabolic risk score for age-related diseases. bioRxiv. https://doi.org/10.1101/2024.11.11.622920.
