# tested on R 3.4.3
# setup required libraries and constants
library(scater)  # tested on version 1.6.3,    install from Bioconductor: source("https://bioconductor.org/biocLite.R"); biocLite("scater")
library(igraph)  # tested on version 1.2.1,   install from CRAN: install.packages("igraph")
library(xgboost) # tested on version 0.6.4.1, install from CRAN: install.packages("xgboost")
library(pROC)    # tested on version 1.12.1,  install from CRAN: install.packages("pROC")
library(caret)   # tested on version 6.0.80,  install from CRAN: install.packages("caret")
library(e1071)   # tested on version 1.6.8,  install from CRAN: install.packages("e1071")
library(xtable)  # tested on version 1.8.2,  install from CRAN: install.packages("xtable")

runCastlePerCellType = function(source, target, name) {
  start_time = Sys.time()
  
  BREAKS=c(-1, 0, 1, 6, Inf)
  nFeatures = 100
 
  
# 1. Load datasets in scater format: loaded files expected to contain "Large SingleCellExperiment" object
  if(class(source)[1]=="character") { source = readRDS(source) }
  if(class(target)[1]=="character") { target = readRDS(target) }
if(class(source)[1]!="SingleCellExperiment") { return }
if(class(target)[1]!="SingleCellExperiment") { return }
ds1 = t(exprs(source)) 
ds2 = t(exprs(target)) 
colnames(ds1) = elementMetadata(source)$feature_symbol
colnames(ds2) = elementMetadata(target)$feature_symbol
sourceCellTypes = as.factor(colData(source)$cell_type1)
targetCellTypes = as.factor(colData(target)$cell_type1)


# remove unlabeled from source
unknownLabels = levels(sourceCellTypes)[grep("not applicable|unclassified|contaminated|unknown", levels(sourceCellTypes))]
if (length(unknownLabels)>0) {
  hasKnownLabel = is.na(match(sourceCellTypes, unknownLabels))
  sourceCellTypes = sourceCellTypes[hasKnownLabel]
  sourceCellTypes = as.factor(as.character(sourceCellTypes))
  ds1 = ds1[hasKnownLabel,]
}

# 2. Unify sets, excluding low expressed genes
source_n_cells_counts = apply(exprs(source), 1, function(x) { sum(x > 0) } )
target_n_cells_counts = apply(exprs(target), 1, function(x) { sum(x > 0) } )
common_genes = intersect( colnames(ds1)[source_n_cells_counts>10], 
                          colnames(ds2)[target_n_cells_counts>10]
)
remove(source_n_cells_counts, target_n_cells_counts)
ds1 = ds1[, colnames(ds1) %in% common_genes]
ds2 = ds2[, colnames(ds2) %in% common_genes]
ds = rbind(ds1[,common_genes], ds2[,common_genes])
isSource = c(rep(TRUE,nrow(ds1)), rep(FALSE,nrow(ds2)))
remove(ds1, ds2)

# 3. Highest mean in both source and target
topFeaturesAvg = colnames(ds)[order(apply(ds, 2, mean), decreasing = T)]

cat("<h2>Source dataset: class distribution</h2>")
print(xtable(t(cbind(table(sourceCellTypes[isSource]),paste0(formatC(prop.table(table(sourceCellTypes[isSource]))*100),"%")))), type="html")
cat("<h2>Target dataset: class distribution</h2>")
print(xtable(t(cbind(table(targetCellTypes),paste0(formatC(prop.table(table(targetCellTypes))*100),"%")))), type="html")

L = length(levels(sourceCellTypes))
summary=data.frame(InstancesInSource=character(L),InstancesInTarget=character(L),Accuracy=character(L),Sensitivity=character(L),Specificity=character(L),AUC=character(L), stringsAsFactors = FALSE)
rownames(summary) = levels(sourceCellTypes)

# for each cell - what is the most probable classification?
targetClassification = as.data.frame(matrix(rep(0,L*sum(!isSource)), nrow=L), row.names = levels(sourceCellTypes))


# iterate over all source cell types
for (cellType in levels(sourceCellTypes)) {
  
  inSourceCellType = as.factor(ifelse(sourceCellTypes == cellType, cellType, paste0("NOT",cellType)))
  
  # 4. Highest mutual information in source
  topFeaturesMi = names(sort(apply(ds[isSource,],2,function(x) { compare(cut(x,breaks=BREAKS),inSourceCellType,method = "nmi") }), decreasing = T))
  
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
  ds0 = model.matrix(~ . , as.data.frame(dsBins[,nUniq>1]))
  remove(dsBins, nUniq)

  cat(paste0("<h2>Classifier for ",cellType,"</h2>"))
  
  inTypeSource = sourceCellTypes == cellType
  # 9. Classify
  xg=xgboost(data=ds0[isSource,] , 
             label=inTypeSource,
             objective="binary:logistic", 
             eta=0.7 , nthread=1, nround=20, verbose=0,
             gamma=0.001, max_depth=5, min_child_weight=10)
             
  # 10. Predict
  inTypeProb = predict(xg, ds0[!isSource, ])
  inTypePred = factor(ifelse(inTypeProb > 0.5,1,0), 0:1)
  
  targetClassification[cellType,] = inTypeProb
  
  # validate against ground truth
  actualInType = factor(ifelse(targetCellTypes == cellType, 1, 0), 0:1) 
  inTypeTargetPrec = mean(targetCellTypes == cellType)
  print(xtable(table(inTypePred, targetCellTypes)), type="html")
  summary[cellType,"InstancesInSource"] = paste0(formatC(mean(inTypeSource)*100 , digits=4),"%")
  summary[cellType,"InstancesInTarget"] = paste0(formatC(inTypeTargetPrec*100 , digits=4),"%")
  cm = confusionMatrix(inTypePred, actualInType)
  summary[cellType,"Accuracy"] = paste0(formatC(as.numeric(cm$overall["Accuracy"])*100 , digits=4),"%")
  summary[cellType,"Sensitivity"] = paste0(formatC(as.numeric(cm$byClass["Sensitivity"])*100 , digits=4),"%")
  
  if (inTypeTargetPrec>0) {
    summary[cellType,"Specificity"] = paste0(formatC(as.numeric(cm$byClass["Specificity"])*100 , digits=4),"%")
    roc_test = roc(as.numeric(actualInType), inTypeProb)
    #plot(roc_test ) 
    summary[cellType,"AUC"] = paste0(formatC(as.numeric(auc(roc_test ))*100 , digits=4),"%")
  } else {
    summary[cellType,"InstancesInTarget"] = "None"
    summary[cellType,"Specificity"] = "-"
    summary[cellType,"AUC"] = "-"
  }
  
}

total_time=as.numeric(difftime(Sys.time(), start_time, units ="secs"))
cat(paste0("Total time: ", round(total_time), " seconds"))

# save summary data
write.csv(cbind(format(t(targetClassification), digits=3),as.character(targetCellTypes)), file = paste0(name,".csv"), row.names = FALSE, quote = FALSE)

# export summary to html
cat("<h2>Summary</h2>")

# Additional questions:
# Of the cell types that were not present in the source datasets, 
#  how many were identified as a specific (wrong) cell type?
hasSomeClassification = apply(targetClassification, 2, function(x) { sum(x>0.5)>0 })
cat(paste0(
  "<p>Target dataset contains ",
  sum(is.na(match(targetCellTypes , levels(sourceCellTypes)))),
  " cells of types that do not exist in source (",
  paste0(levels(targetCellTypes)[is.na(match(levels(targetCellTypes) , levels(sourceCellTypes)))] , collapse = ", "),
  "), of which ",
  sum(hasSomeClassification[is.na(match(targetCellTypes , levels(sourceCellTypes)))]),
  " were classified.</p>"))


# And, of the cell types that were present in the source datasets, 
#  how many were not identified as the same specific cell type? 
# Or of any type?
finalClassification = apply(targetClassification, 2, function(x) { names(which.max(x)) })
hasNoClassification = apply(targetClassification, 2, function(x) { sum(x>0.5)==0 })
finalClassification[hasNoClassification] = "Not Classified"
cat(paste0(
  "<p>Target dataset contains ",
  sum(!is.na(match(targetCellTypes , levels(sourceCellTypes)))),
  " cells of types that exist in source (",
  paste0(levels(targetCellTypes)[!is.na(match(levels(targetCellTypes) , levels(sourceCellTypes)))] , collapse = ", "),
  "), of which ",
  sum(finalClassification[!is.na(match(targetCellTypes , levels(sourceCellTypes)))] != "Not Classified"),
  " were classified, and ",
  sum(finalClassification[!is.na(match(targetCellTypes , levels(sourceCellTypes)))] == targetCellTypes[!is.na(match(targetCellTypes , levels(sourceCellTypes)))]),
  " were classified correctly (",
  paste0(formatC(sum(finalClassification[!is.na(match(targetCellTypes , levels(sourceCellTypes)))] == targetCellTypes[!is.na(match(targetCellTypes , levels(sourceCellTypes)))]) / sum(!is.na(match(targetCellTypes , levels(sourceCellTypes)))) * 100),"%"),
  ").</p>"
))

print(xtable(summary), type="html")


cat("<h2>Cells that got more than one classification</h2>")
hasMoreThanOneClassification = apply(targetClassification, 2, function(x) { sum(x>0.5)>1 })
tmp = t(cbind(
  table(targetCellTypes)[names(table(targetCellTypes[hasMoreThanOneClassification]))],
  table(targetCellTypes[hasMoreThanOneClassification]),
  paste0(formatC(prop.table(table(targetCellTypes[hasMoreThanOneClassification]))*100),"%"),
  paste0(formatC(table(targetCellTypes[hasMoreThanOneClassification])/table(targetCellTypes)[names(table(targetCellTypes[hasMoreThanOneClassification]))]*100),"%")
))
rownames(tmp) = c("Total in cell type","Cells with more than one classification","Out of multi classified","Out of cell type")
print(xtable(tmp), type="html")

cat("<h2>Cells that were not classified</h2>")
hasNoClassification = apply(targetClassification, 2, function(x) { sum(x>0.5)==0 })
tmp = t(cbind(
  table(targetCellTypes)[names(table(targetCellTypes[hasNoClassification]))],
  table(targetCellTypes[hasNoClassification]),
  paste0(formatC(prop.table(table(targetCellTypes[hasNoClassification]))*100),"%"),
  paste0(formatC(table(targetCellTypes[hasNoClassification])/table(targetCellTypes)[names(table(targetCellTypes[hasNoClassification]))]*100),"%")
))
rownames(tmp) = c("Total in cell type","Cells with no classification","Out of no classified","Out of cell type")
print(xtable(tmp), type="html")

}


sink(file="CasTLeClassificationResultsAllTimings.html", append=F)

d1=readRDS("dataset1.rds")
d2=readRDS("dataset2.rds")
runCastlePerCellType(d1, d2, "1_2")
runCastlePerCellType(d2, d1, "2_1")

sink()
