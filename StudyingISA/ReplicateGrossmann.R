# ReplicateGrossman.R: Comparing module stability across 20 runs with same parameters, with imaging features and pathways###
# OUTPUT: Heatmap to view module preservation of rows (pathways) and columns (imaging features)
# Heatmaps of 
# View "track_col_strands" and "track_row_strands" to see how imaging features in modules and pathways in modules
# are preserved, respectively

lapply(list("isa2","eisa","biclust"), require, character.only = TRUE)
dataset = "Data1"
input_dataset <- dataset
import_path <- paste0("~/Documents/MII_Rotation/GrossmannReplication/LeahPipeline/",input_dataset,"_NES_matrix.txt")
nes_data <- as.matrix(read.table(import_path,row.names=1,header=TRUE))
my_normalize <- isa.normalize(nes_data)
####
####



####
####
set.seed(987654321)
unique_collection <- 0
#unique_rows <- 0
#unique_cols <- 0
nruns <- 20
nseeds <- 100
my_direction <- "updown"
# run 10 times
for(i in 1:nruns){

  #isa.option(cor.limit=0.3)
  my_isa <- isa(nes_data,direction=c(my_direction,my_direction),no.seeds = nseeds)
  #every column corresponds to a bicluster
  print(ncol(my_isa$rows))
  my_isa_unique <- isa.unique(my_normalize,my_isa,cor.limit=0.3)
  print(ncol(my_isa_unique$rows))
  #my_biclust <- isa.biclust(my_isa_unique)
  if(i == 1){
    unique_collection <- my_isa_unique
  }
  else {
    unique_collection <- cbind2(unique_collection, my_isa_unique)
  }
 
}

# for grid search!!!
#
#
#
row_strands <- list()
for(module in 1:ncol(unique_collection[1,1][[1]])){
  row_strands[[module]] <- c(1:511)[unique_collection[1,1][[1]][,module] != 0] 
}
col_strands <- list()
for(module in 1:ncol(unique_collection[2,1][[1]])){
  col_strands[[module]] <- c(1:636)[unique_collection[2,1][[1]][,module] != 0] 
}

# here, length colstrands - length rowstrands because we are using all modules in run1 for comparison initially
track_row_strands <- matrix(0,length(col_strands),ncol(unique_collection))
row.names(track_row_strands) <- paste(1,c(1:length(col_strands)),sep=":")
track_col_strands <- matrix(0,length(col_strands),ncol(unique_collection))
row.names(track_col_strands) <- paste(1,c(1:length(col_strands)),sep=":")


#for(compare_run in 4:5){
for(compare_run in 1:ncol(unique_collection)){
  print(paste("run",compare_run))
  #for(compare_module in 1:3){
  for(compare_module in 1:ncol(unique_collection[1,compare_run][[1]])){
    print(paste0("run ",compare_run,"compare_module ",compare_module))
    compare_row <- c(1:511)[unique_collection[1,compare_run][[1]][,compare_module] != 0]
    compare_col <- c(1:636)[unique_collection[2,compare_run][[1]][,compare_module] != 0]
    get_highest_row <- rep(NA,length(row_strands))
    get_highest_col <- rep(NA,length(col_strands))
    
    for(r in 1:nrow(track_row_strands)){
      row_intersect <- intersect(compare_row,row_strands[[r]])
      #get_highest[r] <- (length(row_intersect) + length(col_intersect))/(length(compare_row) + length(compare_col))
      get_highest_row[r] <- (length(row_intersect))/(length(compare_row))
    }

    # is it higher than any max overlap module that came before?
    # 
    if(max(get_highest_row) >= 0.5){
      #print(TRUE)
      for(maxes in which(get_highest_row == max(get_highest_row))){ 
        if(get_highest_row[maxes]>=track_row_strands[maxes,compare_run]){
          track_row_strands[maxes,compare_run] <- max(get_highest_row)
          print(paste0("Row strand most similar",row.names(track_row_strands)[maxes],
                       "highest row overlap",max(get_highest_row)))
        }
        
      }
      #track_row_strands[compare_module,compare_run] <- max(get_highest_row)
    }
    if(max(get_highest_row) < 0.5){
      track_row_strands <- rbind(track_row_strands,rep(0,ncol(unique_collection)))
      track_row_strands[nrow(track_row_strands),compare_run] <- 1
      row.names(track_row_strands)[nrow(track_row_strands)] <- paste(compare_run,compare_module,sep=":")
      # add a new module to the comparison set
      row_strands[[length(row_strands)+1]] = compare_row
    }
    
    for(c in 1:nrow(track_col_strands)){
      col_intersect <- intersect(compare_col,col_strands[[c]])
      #get_highest[r] <- (length(row_intersect) + length(col_intersect))/(length(compare_row) + length(compare_col))
      get_highest_col[c] <- (length(col_intersect))/(length(compare_col))
    }

    #get_highest
    # 
    if(max(get_highest_col) >= 0.5){
      #print(TRUE)
      for(maxes in which(get_highest_col == max(get_highest_col))){
        if(get_highest_col[maxes]>=track_col_strands[maxes,compare_run]){
          track_col_strands[maxes,compare_run] <- max(get_highest_col)
          print(paste0("Col strand most similar",row.names(track_col_strands)[maxes],
                       "highest col overlap",max(get_highest_col)))
        }
      }
      #track_col_strands[compare_module,compare_run] <- max(get_highest_col)
    }
    if(max(get_highest_col) < 0.5){
      track_col_strands <- rbind(track_col_strands,rep(0,ncol(unique_collection)))
      track_col_strands[nrow(track_col_strands),compare_run] <- 1
      row.names(track_col_strands)[nrow(track_col_strands)] <- paste(compare_run,compare_module,sep=":")
      # add a new module to the comparison set
      col_strands[[length(col_strands)+1]] = compare_col
    }
  }
}




# Begin pathway stability heatmap
# get lengths of each row strand and col strand
module_row_lengths <- sapply(row_strands,function(x){length(x)})
row_colors <- brewer.pal(8,"Blues")
row_labels <- c("< 5","5<=x<10","10<=x<25","25<=x<50","50<=x<75","75<=x<100","100<=x<125","125<=x<150")
module_row_lengths_colors <- sapply(module_row_lengths,function(x){
  if(as.integer(x/25)>=1){
    return(row_colors[(as.integer(x/25) + 3)])
  }
  if(x<5){
    return(row_colors[1])
  }
  if(5<=x && x<10){
    return(row_colors[2])
  }
  if(10<=x && x<25){
    return(row_colors[3])
  }  
  })
# get lengths of each col strands

require(gplots)
require(RColorBrewer)
pdf("~/Documents/MII_Rotation/ModuleRobustness/RowModuleStability(Pathways).pdf")
input <- track_row_strands
heatmap.2(input,Rowv=FALSE,Colv=FALSE,trace="none",density="none",xlab="Iterations",ylab="Modules",na.color = "grey",
          col=brewer.pal(9,"Reds"),colsep=0:ncol(input), rowsep=0:nrow(input),sepcolor="black",sepwidth=c(0.02,0.02),
          cellnote=round(input,1),notecex=0.6,notecol = "black",RowSideColors = module_row_lengths_colors)
legend("bottomleft",fill= row_colors,legend = row_labels,cex=0.6,title="Number of pathways")
dev.off()
## end pathway stability heatmap

module_col_lengths <- sapply(col_strands,function(x){length(x)})
col_colors <- brewer.pal(9,"Blues")
col_labels <- c("< 5","5<=x<10","10<=x<25","25<=x<50","50<=x<75","75<=x<100","100<=x<125","125<=x<150","150<=x<175")
module_col_lengths_colors <- sapply(module_col_lengths,function(x){
  if(as.integer(x/25)>=1){
    return(col_colors[(as.integer(x/25) + 3)])
  }
  if(x<5){
    return(col_colors[1])
  }
  if(5<=x && x<10){
    return(col_colors[2])
  }
  if(10<=x && x<25){
    return(col_colors[3])
  }  
})
# get lengths of each col strands

pdf("~/Documents/MII_Rotation/ModuleRobustness/ColModuleStability(Features).pdf")
input <- track_col_strands
heatmap.2(input,Rowv=FALSE,Colv=FALSE,trace="none",density="none",xlab="Iterations",ylab="Modules",na.color = "grey",
          col=brewer.pal(9,"Reds"),colsep=0:ncol(input), rowsep=0:nrow(input),sepcolor="black",sepwidth=c(0.02,0.02),
          cellnote=round(input,1),notecex=0.6,notecol = "black",RowSideColors = module_col_lengths_colors)
legend("bottomleft",fill= col_colors,legend = col_labels,cex=0.6,title="Number of features")
dev.off()
# end imaging feature stability heatmap
#VENNNNN
#install.packages("Vennerable", repos="http://R-Forge.R-project.org")
#install.packages("reshape")
require("reshape")
require("Vennerable")
 # create venn of overlapping modules
intersect(row.names(nes_data)[unique_collection[1,1][[1]][,7]!= 0],row.names(nes_data)[unique_collection[1,4][[1]][,2]!= 0])
intersect(row.names(nes_data)[unique_collection[1,1][[1]][,7]!= 0],row.names(nes_data)[unique_collection[1,5][[1]][,8]!= 0])
#intersect(row.names(nes_data)[unique_collection[1,1][[1]][,7]!= 0],row.names(nes_data)[unique_collection[1,5][[1]][,13]!= 0])
Vstem <-list(row.names(nes_data)[unique_collection[1,1][[1]][,7]!= 0],
             row.names(nes_data)[unique_collection[1,4][[1]][,2]!= 0],
             row.names(nes_data)[unique_collection[1,5][[1]][,8]!= 0])
names(Vstem) <- c("Run 1","Run 4", "Run 5")
myVenn <- Venn(Vstem)
pdf("~/Documents/MII_Rotation/ModuleRobustness/PathwayVenn.pdf")
plot(myVenn,doWeights=FALSE,show=list(Faces=FALSE))
dev.off()
## VENN FEATURES
intersect(colnames(nes_data)[unique_collection[2,1][[1]][,7]!= 0],colnames(nes_data)[unique_collection[2,4][[1]][,2]!= 0])
intersect(colnames(nes_data)[unique_collection[2,1][[1]][,7]!= 0],colnames(nes_data)[unique_collection[2,5][[1]][,8]!= 0])
#intersect(colnames(nes_data)[unique_collection[2,1][[1]][,7]!= 0],colnames(nes_data)[unique_collection[2,5][[1]][,13]!= 0])
Vstem <-list(colnames(nes_data)[unique_collection[2,1][[1]][,7]!= 0],
             colnames(nes_data)[unique_collection[2,4][[1]][,2]!= 0],
             colnames(nes_data)[unique_collection[2,5][[1]][,8]!= 0])
names(Vstem) <- c("Run 1","Run 4", "Run 5")
myVenn <- Venn(Vstem)
pdf("~/Documents/MII_Rotation/ModuleRobustness/FeatureVenn.pdf")
plot(myVenn,doWeights=FALSE,show=list(Faces=FALSE))
dev.off()


## Plotting module stability ##
write.csv(track_strands,paste0("~/Documents/MII_Rotation/ModuleRobustness/TrackStrands",my_direction,"DIR",nruns,"runs",nseeds,"seeds.csv"))
pdf(paste0("~/Documents/MII_Rotation/ModuleRobustness/TrackStrands",my_direction,"DIR",nruns,"runs",nseeds,"seeds.pdf"))
matplot(t(track_strands),type="l")
dev.off()
## Plotting module stability ##

?isa.filter.robust

# on average
mean(rowMeans(track_row_strands)[1:11])
mean(rowMeans(track_col_strands)[1:11])


## Determining the most stable modules ###
# Stable defined as > 0.70 of mean(% intersection) across all runs
good_module_indices <- which(rowMeans(track_strands)>.70) 
# resulting in modules 4,11and 16 being the most highly conserved and larger than 1 pathway, 1 feature
pdf(paste0("~/Documents/MII_Rotation/ModuleRobustness/TrackStrandsBest",my_direction,"DIR",nruns,"runs",nseeds,"seeds.pdf"))
matplot(t(track_strands)[,good_module_indices],type="l")
legend("bottomleft",legend=good_module_indices,col=1:3,lty=1)
dev.off()
## Determining the most stable modules ###

## Dive into best modules
best_modules_rows <- row_strands[good_module_indices]
best_modules_cols <- col_strands[good_module_indices]

best_modules_rows_count_table <- list()
best_modules_cols_count_table <- list()
# r is for one of the three stable modules
for(r in 1:length(best_modules_rows)){
  # count tables for one of the X number of stable modules, where we count how many times a particular index (for feature)
  # of pathway appears in the best matching module for each run
  count_row_table <- matrix(NA,length(best_modules_rows[[r]]),ncol(unique_collection))
  # make the row names the actual indices in these 'stable Modules' for ease of data handling
  row.names(count_row_table) <- best_modules_rows[[r]]
  count_col_table <- matrix(NA,length(best_modules_cols[[r]]),ncol(unique_collection))
  row.names(count_col_table) <- best_modules_cols[[r]]
  
  for(compare_run in 1:ncol(unique_collection)){
    get_highest = rep(NA,ncol(unique_collection[1,compare_run][[1]]))
    for(compare_module in 1:ncol(unique_collection[1,compare_run][[1]])){
      compare_row <- c(1:511)[unique_collection[1,compare_run][[1]][,compare_module] != 0]
      compare_col <- c(1:636)[unique_collection[2,compare_run][[1]][,compare_module] != 0]
      row_intersect <- intersect(compare_row,best_modules_rows[[r]])
      col_intersect <- intersect(compare_col,best_modules_cols[[r]])
      get_highest[compare_module] <- (length(row_intersect) + length(col_intersect))/(length(compare_row) + length(compare_col))
    }
    print(max(get_highest))
    module_matched <- which.max(get_highest)
    count_row_table[,compare_run] <- as.integer(best_modules_rows[[r]] %in% c(1:511)[unique_collection[1,compare_run][[1]][,module_matched] != 0])
    count_col_table[,compare_run] <- as.integer(best_modules_cols[[r]] %in% c(1:636)[unique_collection[2,compare_run][[1]][,module_matched] != 0])
    # get best module across runs, and get it's components
    
  }
  best_modules_rows_count_table[[r]] <- count_row_table
  best_modules_cols_count_table[[r]] <- count_col_table
}
   

row.names(nes_data)[best_modules_rows[[2]][rowSums(best_modules_rows_count_table[[2]]) > 10]]

row.names(nes_data)[best_modules_rows[[3]][rowSums(best_modules_rows_count_table[[3]]) > 10]]
# range of number of runs that contain genes in the module strand
module_no <- c(4,11,16)
for(c in c(1:3)){
  print(range(rowSums(best_modules_rows_count_table[[c]])))
  # write down pathways and features with appearance in more than 10/20 runs
  write(paste0(row.names(nes_data)[best_modules_rows[[c]][rowSums(best_modules_rows_count_table[[c]]) > 10]],"\n"),
        paste0("~/Documents/MII_Rotation/ModuleRobustness/HighlyConservedPathFeat",
               my_direction,"DIR",nruns,"runs",nseeds,"seeds",module_no[c],".txt"))
  write("\n Features",
        paste0("~/Documents/MII_Rotation/ModuleRobustness/HighlyConservedPathFeat",
               my_direction,"DIR",nruns,"runs",nseeds,"seeds",module_no[c],".txt"),append=TRUE)
  write(paste0(colnames(nes_data)[best_modules_cols[[c]][rowSums(best_modules_cols_count_table[[c]]) > 10]],"\n"),
        paste0("~/Documents/MII_Rotation/ModuleRobustness/HighlyConservedPathFeat",
               my_direction,"DIR",nruns,"runs",nseeds,"seeds",module_no[c],".txt"),append=TRUE)
  
}


