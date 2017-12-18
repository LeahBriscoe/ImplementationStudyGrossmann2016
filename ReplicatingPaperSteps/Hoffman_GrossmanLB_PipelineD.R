#PIPELINE D: Use Iterative Signature Algorithm to produce pathway-imaging modules
#             Validate modules on Dataset 2 using permutation test and plot valid modules in heatmap
#             Also write to file the features and pathways in these validated modules
#             One Heatmap per Module
#OUTPUT: PDFs of Modules, txt files with moduled components


#args <- commandArgs(trailingOnly = TRUE)
#dataset = args[1]
dataset = "Data1"
# Biclustering 
#source("http://bioconductor.org/biocLite.R")
#biocLite("reactome.db")
#install.packages('isa2')
#install.packages('biclust')
#install.packages('AnnotationDbi')
my_module_sizes
lapply(list("isa2","eisa","biclust","Biobase",'org.Hs.eg.db',"AnnotationDbi","methods","reactome.db"), require, character.only = TRUE)
input_dataset <- dataset
import_path <- paste0("~/Documents/MII_Rotation/ModuleRobustness/",input_dataset,"_NES_matrix.txt")
nes_data <- as.matrix(read.table(import_path,row.names=1,header=TRUE))
nes_data_df <- read.table(import_path,row.names=1,header=TRUE)
input_dataset <- "Data2"
import_path <- paste0("~/Documents/MII_Rotation/ModuleRobustness/",input_dataset,"_NES_matrix.txt")
nes_data2 <- as.matrix(read.table(import_path,row.names=1,header=TRUE))
corrected_d2 <- nes_data2[,colnames(nes_data)]
nes_data2 <- corrected_d2

dir.create("~/Documents/MII_Rotation/ModuleRobustness/ISA_Modules_Results/")
#ncol(nes_data2)
#636 imaging features in common between dataset 1 and dataset 2
#length(intersect(colnames(nes_data2),colnames(nes_data)))
# Give up on expressionset method with ISA for now because of database issue
# my_eset <- ExpressionSet(nes_data,annotation='reactome')
# my_isa_eset <- ISANormalize(my_eset, prenormalize = FALSE)
# my_eisa <- ISA(my_isa_eset,flist=NA,uniqueEntrez = FALSE,thr.gene=seq(1.5,2.5,0.5),thr.cond=seq(1.5,2.5,0.5))
# my_isa <- isa(nes_data,thr.row = seq(1.5,2.5,0.5),thr.col = seq(1.5,2.5,0.5),direction = c('updown',"updown"))
# if (interactive()) {
#   sapply(my_isa, function(b) image(outer(my_isa$rows[,b],
#                                        my_isa$columns[,b]),
#                                  main=paste("Module", b)))  
# }

# returns two matrices, one is a transpose of the others

my_seeds <- 987654321
my_direction <- "down"
my_isa <- isa(nes_data,thr.row=seq(1.5,2.5,0.5),thr.col=seq(1.5,2.5,0.5),no.seeds=my_seeds,direction=c(my_direction,my_direction))
#every column corresponds to a bicluster
print(ncol(my_isa$row))
print(ncol(my_isa$col))
# I should just do this!
my_normalize <- isa.normalize(nes_data)
my_isa_unique <- isa.unique(my_normalize,my_isa,cor.limit=0.3)
print(ncol(my_isa_unique$row))
print(ncol(my_isa_unique$col))
# updown,updown, 100 seeds = 25 modules, 5 modules after unique
# updown, updown, 636 seeds = 63 modules, 15 modules afer unique 
# updown, updown, 511seeds = 51 modules, 11 modules afer unique 
# up, up, 100 seeds = 39, 12
# up up 511 seeds = 68, 13
# up up 636 seeds = 59, 10
# updown,updown,1147 seeds = 18
#show sie of modules
my_module_sizes <- cbind( colSums(my_isa_unique$rows!=0), colSums(my_isa_unique$columns!=0) )
#  correlation of isa and make histogram
cc <- cor(my_isa_unique$rows)
pdf(paste("~/Documents/MII_Rotation/ModuleRobustness/ISA_Modules_Results/",
          ncol(my_isa_unique$row),my_seeds,my_direction,".pdf",sep="_"))
if (interactive()) { hist(cc[lower.tri(cc)],10,main=paste0("Correlation of ", ncol(my_isa_unique$row),
                          " Modules with ",my_seeds," seeds, ",my_direction," direction"),cex.main=1) }
dev.off()
number_modules <- ncol(my_isa_unique$row)
#plotModules(my_isa_unique,data=nes_data)
# Convert to biclust object which binarizes the numbers to true or false in module
my_biclust <- isa.biclust(my_isa_unique)
#my_obiclust <- isa.biclust(my_isa)
###Remove NA from NES matrix
new_nes_data <- apply(nes_data,2,function(x){
  test <- x
  test[is.na(x)] <- 0
  return(test)})
sum(apply(new_nes_data,1,function(x){is.na(x)}))
###
### Dataset1 correlation matrix
d1_cor <- matrix(NA,number_modules,2)
colnames(d1_cor) <- c("C_x_features","C_y_pathways")
d2_cor <- matrix(NA,number_modules,2)
colnames(d2_cor) <- c("C_x_features","C_y_pathways")
for(i in 1:ncol(my_isa_unique$rows)){
  print(i)
  my_module_correlation <- as.matrix(nes_data[(my_biclust@RowxNumber[,i]),(my_biclust@NumberxCol[i,])])
  row.names(my_module_correlation) <- row.names(nes_data)[(my_biclust@RowxNumber[,i])]
  colnames(my_module_correlation) <- colnames(nes_data)[(my_biclust@NumberxCol[i,])]
  #mean of whole correlation matrix
  # mean of correlation of pairs of features in a module
  C_x <- mean(cor(my_module_correlation,use="complete.obs",method="spearman"))
  feature_correlation_means <- rowMeans((cor(my_module_correlation,use="complete.obs",method="spearman")))
  print(which(feature_correlation_means==max(feature_correlation_means)))
  # mean of correlation of pairs of pathways in a module
  C_y <- mean(cor(t(my_module_correlation),use="complete.obs",method="spearman"))
  pathway_correlation_means <- rowMeans((cor(t(my_module_correlation),use="complete.obs",method="spearman")))
  #print(which(pathway_correlation_means==max(pathway_correlation_means)))
  print(names(tail(sort(pathway_correlation_means),3)))
  
  d1_cor[i,] <- c(C_x,C_y)
  d2_test <- as.matrix(nes_data2[row.names(my_module_correlation),colnames(my_module_correlation)])
  C_x <- mean(cor(d2_test,use="complete.obs",method="spearman"))
  # mean of correlation of pairs of pathways in a module
  C_y <- mean(cor(t(d2_test),use="complete.obs",method="spearman"))
  d2_cor[i,] <- c(C_x,C_y)
}
true_d1r <- apply(d1_cor,1,sum)

# Permutation test
collect_d2r <- matrix(NA,1000,number_modules)
for(m in 1:nrow(my_module_sizes)){
  #rows pathways, col features
  # 1st size is for pathways
  for(p in 1:1000){
    pathway_sample <- sample(nrow(nes_data2),my_module_sizes[m,1])
    feature_sample <- sample(ncol(nes_data2),my_module_sizes[m,2])
    d2_sample <- as.matrix(nes_data2[pathway_sample,feature_sample])
    collect_d2r[p,m] <- mean(cor(d2_sample,use="complete.obs",method="spearman")) + 
                        mean(cor(t(d2_sample),use="complete.obs",method="spearman"))
  }
}
pdf(paste("~/Documents/MII_Rotation/ModuleRobustness/ISA_Modules_Results/HistD2",
          ncol(my_isa_unique$row),my_seeds,my_direction,".pdf",sep="_"))
hist(collect_d2r[,4],breaks=100)
dev.off()

permutation_pval <- matrix(NA,number_modules,1)
for(t in 1:length(true_d1r)){
  permutation_pval[t,1] <- sum(collect_d2r[,t]>true_d1r[t])/1000
}
(sum(permutation_pval < 0.05,na.rm=TRUE))

permutation_pval_adjust <- p.adjust(permutation_pval,method="BH")
(sum(permutation_pval_adjust < 0.05,na.rm=TRUE))
valid_modules = NA
for(n in 1:length(permutation_pval_adjust)){
  if(permutation_pval_adjust[n] < 0.05 && !is.na(permutation_pval_adjust[n])){
    valid_modules = c(valid_modules,n)
    print(n)
    write(paste("Number Pathways:",sum(my_biclust@RowxNumber[,n])),file=paste("~/Documents/MII_Rotation/ModuleRobustness/ISA_Modules_Results/Module"
                                                                      ,n,ncol(my_isa_unique$row),my_seeds,my_direction,".txt",sep="_"),append=FALSE)
    write(paste("Number Features:",sum(my_biclust@NumberxCol[n,])),file=paste("~/Documents/MII_Rotation/GrossmannReplication/LeahPipeline/ISA_Modules_Results/Module"
                                                                              ,n,ncol(my_isa_unique$row),my_seeds,my_direction,".txt",sep="_"),append=TRUE)
    
    write(row.names(nes_data)[(my_biclust@RowxNumber[,n])],file=paste("~/Documents/MII_Rotation/ModuleRobustness/ISA_Modules_Results/Module"
                                                                      ,n,ncol(my_isa_unique$row),my_seeds,my_direction,".txt",sep="_"),append=TRUE)
    write("\n",file=paste("~/Documents/MII_Rotation/GrossmannReplication/LeahPipeline/ISA_Modules_Results/Module"
                          ,n,ncol(my_isa_unique$row),my_seeds,my_direction,".txt",sep="_"),append=TRUE)
    write(paste(colnames(nes_data)[(my_biclust@NumberxCol[n,])]),file=paste("~/Documents/MII_Rotation/ModuleRobustness/ISA_Modules_Results/Module"
                                                                            ,n,ncol(my_isa_unique$row),my_seeds,my_direction,".txt",sep="_"),append=TRUE)
  }
}
valid_modules = valid_modules[-1]
library(ggplot2)
library(reshape2)
library(gplots)


for(v in valid_modules){
  if( sum(my_biclust@RowxNumber[,v]) > 1 && sum(my_biclust@NumberxCol[v,]) > 1){
    print(v)
    my_module = nes_data[my_biclust@RowxNumber[,v],my_biclust@NumberxCol[v,]]
    scale_my_module <- scale(t(my_module))
    colnames(scale_my_module) <- gsub("REACTOME_","",colnames(scale_my_module))
    my_module.m = melt(scale_my_module)
    names(my_module.m) <- c("Features","Pathways","NES")
    pdf(paste("~/Documents/MII_Rotation/ModuleRobustness/ISA_Modules_Results/GGPLOT_Module"
              ,v,ncol(my_isa_unique$row),my_seeds,my_direction,".pdf",sep="_"),width=8, height=8)
    my_plot <- ggplot(my_module.m, aes(Pathways, Features,fill=NES)) +
      geom_tile(colour = "white") +
      scale_fill_gradient(low = "green", high = "blue") +
      ylab("List of genes ") +
      xlab("List of patients") +
      theme(text = element_text(size=8),
          # legend.title = element_text(size = 6),
          #   legend.text = element_text(size = 6),
          #   plot.title = element_text(size=8),
          #   axis.title=element_text(size=6,face="bold"),
            axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(fill = "Expression level")
    print(my_plot)
    dev.off()
    my_col <- colorRampPalette(c("green","blue"))(n=100)
    feature_colors <- unlist(lapply(row.names(scale_my_module),function(x){
      if(grepl("LoG",x)) '#FFC0CB' #pink
      else if(grepl('Wavelet_',x)) '#FF0000' # red
      else if(grepl('GLSZM_|GLCM|RLGL',x)) '#EB8922' #orange texture?
      else if(grepl('Shape_',x)) '#7D0FB4'
      else if(grepl('Stats',x)) '#000000' #black
    }))
    pdf(paste("~/Documents/MII_Rotation/ModuleRobustness/ISA_Modules_Results/Heatmap_Module"
              ,v,ncol(my_isa_unique$row),my_seeds,my_direction,".pdf",sep="_"),width=8, height=8)
  
    heatmap.2(scale_my_module,col=my_col,trace="none", density="none",margins=c(25,15),RowSideColors = feature_colors,cexCol = 0.8)
    legend("topright",legend=c("LoG","Wavelet","Texture","Shape","Stats"),fill=c('#FFC0CB','#FF0000','#EB8922','#7D0FB4','#000000'),cex=0.7)
    dev.off()
    #dim(scale_my_module)
    #length(feature_colors)
    
  }
}
# # Assess module overlap
# for(v_1 in valid_modules){
#   for(v_2 in valid_modules){
#     features1 <- colnames(nes_data)[(my_biclust@NumberxCol[v_1,])]
#     #replacer <- c("LoG_","Wavelet_","GLCM_","GLSZM_","Shape_","Stats_"
#     for(r in c("LoG_","Wavelet_","GLCM_","GLSZM_","Shape_","Stats_")){
#       features1 <- gsub(replacer,"",features1)
#     }
#     features2 <- colnames(nes_data)[(my_biclust@NumberxCol[v_2,])]
#     print(paste(v_1,v_2))
#     print(length(intersect(features1,features2)))
#   }
# }
# grepl("Wavelet_[[:alpha:]]{3}_",features1)
# gsub("Wavelet_[[:upper:]]{3}_","",features1)
#drawHeatmap2(new_nes_data,my_obiclust,plotAll = TRUE)
#bubbleplot(new_nes_data,my_biclust,showLabels=TRUE)
#data <- isa.in.silico()
#typeof(data[[1]])


#Attempt to merge modules myself?
# my_normalize <- isa.normalize(nes_data)
# my_row_seeds <- generate.seeds(nrow(nes_data))
# my_col_seeds <- generate.seeds(ncol(nes_data))
# for(i in seq(1.5,2.5,0.5)){
#   for(j in seq(1.5,2.5,0.5)){
#     my_iterate <- isa.iterate(my_normalize,row.seeds=my_row_seeds,col.seeds=my_col_seeds,
#                               thr.row = i,thr.col = j,direction = c('updown',"updown"))
#     my_unique <- isa.unique(my_normalize,my_iterate,cor.limit=0.3)
#     print(paste("before unique",ncol(my_iterate$rows)," after unique: ",ncol(my_unique$rows)))
#     my_filter <- isa.filter.robust(normed.data = my_normalize, data = nes_data,isares = my_unique)
#     print(ncol(my_filter$row))
#     print(ncol(my_filter$col))
#     (my_filter$rundata$corx)
#     
#   }
# }

# make one object
# working with module 1 here
# 1.5 to 2.5 by 0.5 for row and coilumn thresholds
# workflow: isa.normalize, generate.seeds, isa.iterate, isa.unique, isa.filter, robust,merging similar modules

