# PIPELINE-B: For each column of weighted correlation matrix from PIPELINE A, corresponding to 1 imaging feature
#             create ranked gene list based on weighted correlation values: 1 column of gene id, 1 column of weights
# OUTPUT: Numerous ranked gene lists: one per imaging feature across all genes

dataset="Data1" #or Data2 for dataset 2
#lapply(list("biomaRt","annotate","org.Hs.eg.db"), require, character.only = TRUE)
#spearman_weighted <- read.csv(paste0("/u/home/b/briscoel/scratch/GrossmanLB_Pipeline/",dataset,"save_weighted_spearman.csv"),header=TRUE, row.names=1)
spearman_weighted <- read.csv(paste0("~/Documents/MII_Rotation/ModuleRobustness/UnofficialData/",dataset,"save_weighted_spearman.csv"),header=TRUE, row.names=1)

#export_path <- paste0("/u/home/b/briscoel/scratch/GrossmanLB_Pipeline/",dataset,"EntrezGeneWeightsPerRadiomicFeature/")
export_path <- paste0("~/Documents/MII_Rotation/ModuleRobustness/UnofficialData/",dataset,"EntrezGeneWeightsPerRadiomicFeature/")

dir.create(export_path)

# substitute "geneid." in row names with ""
genenames_before_conversion <- gsub("geneid.","",row.names(spearman_weighted))
# convert to gene symbol
for(n in 1:ncol(spearman_weighted)){
  safe_list1 <- as.matrix(spearman_weighted[!is.infinite(spearman_weighted[,n]),n])
  
  
  safe_list2 <- cbind(genenames_before_conversion[!is.infinite(spearman_weighted[,n])], safe_list1) 
  print("NA" %in% genenames_before_conversion[!is.infinite(spearman_weighted[,n])])
  safe_list3 <- rbind(c("GeneName","Weight"),safe_list2)
  
  write.table(safe_list3,paste(export_path,colnames(spearman_weighted)[n],".rnk",sep=""),
              sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
  
}
