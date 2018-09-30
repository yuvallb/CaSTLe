# CaSTLe
## Classification of Single cells by Transfer Learning

Label cells in an Single cell RNA sequencing (scRNA-seq) experiment, using knowledge learnt from other experiments on similar cell types.
R code for loading the Source and Target datasets, feature selection, transformation and classification using XGBoost.
  
It is assumed that the Source and Target datasets are available in files in scater format of “Large SingleCellExperiment” object, regarded as sourceDataset.rds and tragetDataset.rds in the code.

The R code uses the following libraries:
  * scater: Single-cell analysis toolkit for gene expression data in R (McCarthy et al., 2017) is used for loading and reading scRNA-seq datasets, in a well-defined format.
  * Xgboost: eXtreme Gradient Boosting (Chen et al. , 2017) is used for building the classification model itself.
  * Igraph (Csardi & Nepusz, 2006) is used for calculating mutual information.

