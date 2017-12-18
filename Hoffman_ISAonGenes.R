# Run Iterative Signature Algorithm on Gene Expressino data directly 

install.packages('isa2',dependencies=TRUE,repos='http://cran.r-project.org')
lapply(list("isa2","biclust"), require, character.only = TRUE)
dataset <- "CleanData1"
#import_path <- paste0("/u/home/b/briscoel/scratch/GrossmanLB_Pipeline/ImplementationStudy/",dataset,"Expression.csv")
import_path <- paste0("~/Documents/MII_Rotation/ModuleRobustness/UnofficialData/",dataset,"Expression.csv")
#export_path <- paste0("/u/home/b/briscoel/scratch/GrossmanLB_Pipeline/ImplementationStudy/",dataset)
export_path <- paste0("~/Documents/MII_Rotation/ModuleRobustness/",dataset)

edat <- as.matrix(read.csv(import_path,row.names=1,header=TRUE,check.names=FALSE))

set.seed(987654321)
nseeds <- 100
my_direction <- "updown"

my_gene_isa <- isa(edat,direction=c(my_direction,my_direction),no.seeds = nseeds)
module_columns <- my_gene_isa$columns
module_rows <- my_gene_isa$rows
write.csv(module_columns,paste0(export_path,"ModuleColumns.csv"))
write.csv(module_rows,paste0(export_path,"ModuleRows.csv"))

