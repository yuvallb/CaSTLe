# setup required libraries and constants
library(scater)
library(xgboost)
library(igraph)
BREAKS=c(-1, 0, 1, 6, Inf)
nFeatures = 100

# 1. Load datasets in scater format: loaded files expected to contain "Large SCESet" object
source = readRDS("sourceDataset.rds")
target = readRDS("tragetDataset.rds")
ds1 = t(exprs(source)) 
ds2 = t(exprs(target)) 
labels = as.factor(pData(source)[,"cell_type1"])

# 2. Unify sets, excluding low expressed genes
common_genes = intersect( fData(source)[fData(source)[,"n_cells_exprs"]>10,"feature_symbol"], 
fData(target)[fData(target)[,"n_cells_exprs"]>10,"feature_symbol"]
)
ds1 = ds1[, colnames(ds1) %in% common_genes]
ds2 = ds2[, colnames(ds2) %in% common_genes]
ds = rbind(ds1[,common_genes], ds2[,common_genes])
isSource = c(rep(TRUE,nrow(ds1)), rep(FALSE,nrow(ds2)))
remove(ds1, ds2)

# 3. Highest mean in both source and target
topFeaturesAvg = colnames(ds)[order(apply(ds, 2, mean), decreasing = T)]

# 4. Highest mutual information in source
topFeaturesMi = names(sort(apply(ds[isSource,],2,function(x) { compare(cut(x,breaks=BREAKS),labels,method = "nmi") }), decreasing = T))

# 5. Top n genes that appear in both mi and avg
selectedFeatures = union(head(topFeaturesAvg, nFeatures) , head(topFeaturesMi, nFeatures) )

# 6. remove correlated features
tmp = cor(ds[,selectedFeatures], method = "pearson")
tmp[!lower.tri(tmp)] = 0
selectedFeatures = selectedFeatures[apply(tmp,2,function(x) any(x < 0.9))]
remove(tmp)

# 7,8. Convert data from continous to binned dummy vars
# break datasets to bins
dsBins = apply(ds[, selectedFeatures], 2, cut, breaks= BREAKS)
# use only bins with more than one value
nUniq = apply(dsBins, 2, function(x) { length(unique(x)) })
# convert to dummy vars
ds = model.matrix(~ . , as.data.frame(dsBins[,nUniq>1]))
remove(dsBins, nUniq)

# 9. Classify
train = runif(nrow(ds[isSource,]))<0.8
# slightly different setup for multiclass and binary classification
if (length(unique(labels)) > 2) {
  xg=xgboost(data=ds[isSource,][train, ] , 
       label=as.numeric(labels[train])-1,
       objective="multi:softmax", num_class=length(unique(labels)),
       eta=0.7 , nthread=5, nround=20,
       gamma=0.001, max_depth=5, min_child_weight=10)
} else {
  xg=xgboost(data=ds[isSource,][train, ] , 
       label=as.numeric(labels[train])-1,
       eta=0.7 , nthread=5, nround=20,
       gamma=0.001, max_depth=5, min_child_weight=10)
}

# 10. Predict
predictedClasses = predict(xg, ds[!isSource, ])


