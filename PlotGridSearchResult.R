require(gplots)
require(RColorBrewer)
setwd("~/Documents/MII_Rotation/GrossmannReplication/ImplementationStudy/")
file_names <- list.files(path=".",pattern="NES_.",full.names=F)
clean_file_names <- gsub(".csv","",file_names)
for(f in 1:length(file_names)){
  my_matrix <- read.csv(file_names[f],row.names=1)
  row.names(my_matrix) <- seq(0.5,5,0.5) 
  colnames(my_matrix) <- seq(100,5000,200) 
  pdf(paste0("NES-plots/",clean_file_names[f],".pdf"))
  input <- t(as.matrix(my_matrix))
  heatmap.2(input,Rowv=FALSE,Colv=FALSE,trace="none",density="none",xlab="ISA Threshold",ylab="No. Seeds",na.color = "grey",
            col=brewer.pal(9,"Reds"),colsep=0:ncol(input), rowsep=0:nrow(input),sepcolor="black",sepwidth=c(0.02,0.02),
            cellnote=round(input,1),notecex=0.6,notecol = "black")
  dev.off()
}







