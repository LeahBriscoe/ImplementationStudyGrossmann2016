# Hoffman_ISAonGenes: Pre-GSEA data prep, unfinished code. I intended to use resulting weights from running ISA on genes
#                     and then looking at enrichment of different pathways in the modules that form
# Output gene set for each patient in dataset 1 for use in GSEA
require(gplots)
require(RColorBrewer)
setwd("~/Documents/MII_Rotation/ModuleRobustness/UnofficialData/")
file_names <- list.files(path=".",pattern="CleanData1Module",full.names=F)
my_rows <- read.csv(file_names[2],row.names=1)
my_cols <- read.csv(file_names[1],row.names=1)
#running the whole pipeline on genes produces 292 modules - output each gene list
edata <- read.csv("~/Documents/MII_Rotation/ModuleRobustness/UnofficialData/CleanData1Expression.csv",row.names=1)
gene_names <- row.names(edata)

for(n in 1:ncol(spearman_weighted)){
  safe_list1 <- as.matrix(spearman_weighted[!is.infinite(spearman_weighted[,n]),n])
  
  
  safe_list2 <- cbind(genenames_before_conversion[!is.infinite(spearman_weighted[,n])], safe_list1) 
  print("NA" %in% genenames_before_conversion[!is.infinite(spearman_weighted[,n])])
  safe_list3 <- rbind(c("GeneName","Weight"),safe_list2)
  
  write.table(safe_list3,paste(export_path,colnames(spearman_weighted)[n],".rnk",sep=""),
              sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
  
}