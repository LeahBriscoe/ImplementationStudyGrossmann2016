# Test ISA at multiple threshold levels and seed numbers
# Output: CSV Matrix with number of modules developed at each parameter set (various threshold levels and seed numbers)
#         - Modules after isa.iterate
#         - Modules after isa.unique
#         - Modules after isa.filter.robust


install.packages('isa2',dependencies=TRUE,repos='http://cran.r-project.org')
lapply(list("isa2","biclust"), require, character.only = TRUE)
dataset = "Data1"
input_dataset <- dataset
import_path <- paste0("~/Documents/MII_Rotation/ModuleRobustness/",input_dataset,"_NES_matrix.txt")
nes_data <- as.matrix(read.table(import_path,row.names=1,header=TRUE))
my_normalize <- isa.normalize(nes_data)


tune_thresholds <- seq(0.5,5,0.5) 
tune_seeds <- seq(100,5000,200) 

collect_nmodules_isa <- matrix(NA,length(tune_thresholds),length(tune_seeds))
collect_nmodules_unique <- matrix(NA,length(tune_thresholds),length(tune_seeds))
collect_nmodules_filter <- matrix(NA,length(tune_thresholds),length(tune_seeds))
collect_avgrowsizemodules_filter <- matrix(NA,length(tune_thresholds),length(tune_seeds))
collect_avgcolsizemodules_filter <- matrix(NA,length(tune_thresholds),length(tune_seeds))

for(t_s in 1:length(tune_seeds)){
  my_row_seeds <- generate.seeds(nrow(nes_data),count=tune_seeds[t_s])
  my_col_seeds <- generate.seeds(ncol(nes_data),count=tune_seeds[t_s])
  for(t_t in 1:length(tune_thresholds)){
    my_iterate <- isa.iterate(my_normalize,row.seeds=my_row_seeds,col.seeds=my_col_seeds,
                              thr.row = tune_thresholds[t_t],thr.col = tune_thresholds[t_t],direction = c('updown',"updown"))
    my_unique <- isa.unique(my_normalize,my_iterate,cor.limit=0.3)
    my_filter <- isa.filter.robust(normed.data = my_normalize, data = nes_data,isares = my_unique)
    
    collect_nmodules_isa[t_t,t_s] <- ncol(my_iterate$rows)
    collect_nmodules_unique[t_t,t_s] <- ncol(my_unique$rows)
    collect_nmodules_filter[t_t,t_s] <- ncol(my_filter$rows)
    binarized_r <- ceiling(abs(my_filter$rows))
    collect_avgrowsizemodules_filter[t_t,t_s] <- mean(colSums(binarized_r)) 
    binarized_c <- ceiling(abs(my_filter$columns))
    collect_avgcolsizemodules_filter[t_t,t_s] <- mean(colSums(binarized_c)) 

  }
}

write.csv(collect_nmodules_isa,"~/Documents/MII_Rotation/ModuleRobustness/NES_NumberModules_ISAiterate.csv")
write.csv(collect_nmodules_unique,"~/Documents/MII_Rotation/ModuleRobustness/NES_NumberModules_ISAunique.csv")
write.csv(collect_nmodules_filter,"~/Documents/MII_Rotation/ModuleRobustness/NES_NumberModules_ISAfilter.csv")
write.csv(collect_avgrowsizemodules_filter,"~/Documents/MII_Rotation/ModuleRobustness/NES_SizeRowModules_ISAfilter.csv")
write.csv(collect_avgcolsizemodules_filter,"~/Documents/MII_Rotation/ModuleRobustness/NES_SizeColModules_ISAfilter.csv")