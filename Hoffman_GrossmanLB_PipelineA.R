dataset <- "Data1" # or "Data2" for dataset2

export_path <- "~/Documents/MII_Rotation/ModuleRobustness/UnofficialData/"
#export_path <- "/u/home/b/briscoel/scratch/GrossmanLB_Pipeline/"
idat <- read.csv(paste0(export_path,dataset,"Imaging.csv"),row.names=1,header=TRUE,check.names=FALSE)
edat <- read.csv(paste0(export_path,dataset,"Expression.csv"),row.names=1,header=TRUE,check.names=FALSE)
#cdat <- read.csv(paste0(export_path,dataset,"Clinical.csv"),row.names=1,header=TRUE,check.names=FALSE)
# sort all
sort_idat <- idat[order(rownames(idat)),]
sort_edat <- edat[,order(colnames(edat))]
#sort_cdat <- cdat[order(rownames(cdat)),]
idat <- sort_idat[,4:ncol(idat)]
edat <- sort_edat
#cdat <- sort_cdat

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
spearman_weighted <- (spearman_rho) * (spearman_pvalue)
row.names(spearman_weighted) <- row.names(edat)
spearman_weighted <- spearman_weighted[,-1]
write.csv(spearman_weighted,paste0(export_path,dataset,"save_weighted_spearman.csv"),quote=FALSE)

# save_weighted_spearman <- save_spearman * save_spearman_pvalue

#write.csv(save_spearman_pvalue,"/u/home/b/briscoel/scratch/save_spearman_pvalue.csv",quote=FALSE)
#write.csv(save_weighted_spearman,"/u/home/b/briscoel/scratch/save_weighted_spearman.csv",quote=FALSE)