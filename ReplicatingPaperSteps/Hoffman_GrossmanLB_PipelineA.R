# PIPELINE-A: Calculated weighted correlation values for imaging and gene expression data 
# OUTPUT: Correlation matrix and Weighted correlation matrix to UnofficialData folder

# Dataset 1 is the train dataset and dataset is the test dataset for the methods in this paper
dataset <- "Data1" # or "Data2" for dataset2

export_path <- "~/Documents/MII_Rotation/ModuleRobustness/UnofficialData/"
# Hoffman2 path:
#export_path <- "/u/home/b/briscoel/scratch/GrossmanLB_Pipeline/"
idat <- read.csv(paste0(export_path,dataset,"Imaging.csv"),row.names=1,header=TRUE,check.names=FALSE)
edat <- read.csv(paste0(export_path,dataset,"Expression.csv"),row.names=1,header=TRUE,check.names=FALSE)
#cdat <- read.csv(paste0(export_path,dataset,"Clinical.csv"),row.names=1,header=TRUE,check.names=FALSE)
# sort all
sort_idat <- idat[order(rownames(idat)),]
sort_edat <- edat[,order(colnames(edat))]
idat <- sort_idat[,4:ncol(idat)] # remove unwanted columns
edat <- sort_edat


save_spearman <- apply(idat,2,function(x)
{apply(edat,1,function(y) {core_result <- cor.test(as.numeric(x),as.numeric(y),method="spearman")
  core_result1 <- core_result$estimate
  core_result2 <- -log10(core_result$p.value)
  c(core_result1,core_result2)
  })
})
write.csv(save_spearman,paste0(export_path,dataset, "save_spearman.csv"),quote=FALSE)
X <- c(1:nrow(save_spearman))
save_spearman <- cbind(X,save_spearman)
spearman_rho <- save_spearman[(save_spearman[,1] %% 2 == 1),]
spearman_pvalue <- save_spearman[(save_spearman[,1] %% 2 == 0),]
# multiply negative log 10(pvalue) by correlation value
spearman_weighted <- (spearman_rho) * (spearman_pvalue)
row.names(spearman_weighted) <- row.names(edat)
spearman_weighted <- spearman_weighted[,-1]
write.csv(spearman_weighted,paste0(export_path,dataset,"save_weighted_spearman.csv"),quote=FALSE)
