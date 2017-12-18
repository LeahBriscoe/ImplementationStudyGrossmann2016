# ISAonGenes: Look at module preservation from two runs of Hoffman_ISAonGenes.R
# No file exported, view "track_col_strands" and "track_row_strands" to see how patients in modules and genes in modules
# are preserved, respectively
require(gplots)
require(RColorBrewer)
setwd("~/Documents/MII_Rotation/GrossmannReplication/ImplementationStudy/")
file_names_run1 <- list.files(path="./First_ISAonGenes/",pattern="CleanData1Module",full.names=T)
my_rows_1 <- read.csv(file_names_run1[2],row.names=1)
my_cols_1 <- read.csv(file_names_run1[1],row.names=1)
file_names_run2 <- list.files(path="./Second_ISAonGenes/",pattern="CleanData1Module",full.names=T)
my_rows_2 <- read.csv(file_names_run2[2],row.names=1)
my_cols_2 <- read.csv(file_names_run2[1],row.names=1)
#running the whole pipeline on genes produces 292 modules
unique_collection <-list()
unique_collection[[1]] <- list(my_rows_1,my_cols_1)
unique_collection[[2]] <- list(my_rows_2,my_cols_2)

row_strands <- list()
for(module in 1:ncol(unique_collection[[1]][[1]])){
  row_strands[[module]] <- c(1:21766)[unique_collection[[1]][[1]][,module] != 0] 
}
col_strands <- list()
for(module in 1:ncol(unique_collection[[1]][[2]])){
  col_strands[[module]] <- c(1:262)[unique_collection[[1]][[2]][,module] != 0] 
}


# here, length colstrands - length rowstrands because we are using all modules in run1 for comparison initially
track_row_strands <- matrix(0,length(col_strands),length(unique_collection))
row.names(track_row_strands) <- paste(1,c(1:length(col_strands)),sep=":")
track_col_strands <- matrix(0,length(col_strands),length(unique_collection))
row.names(track_col_strands) <- paste(1,c(1:length(col_strands)),sep=":")

c(1:21766)[unique_collection[[1]][[1]][,5] != 0]
c(1:21766)[unique_collection[[2]][[1]][,5] != 0]

#for(compare_run in 4:5){
for(compare_run in 1:length(unique_collection)){
  print(paste("run",compare_run))
  #for(compare_module in 1:3){
  for(compare_module in 1:ncol(unique_collection[[compare_run]][[1]])){
    print(paste0("run ",compare_run,"compare_module ",compare_module))
    compare_row <- c(1:21766)[unique_collection[[compare_run]][[1]][,compare_module] != 0] 
    compare_col <- c(1:262)[unique_collection[[compare_run]][[2]][,compare_module] != 0] 
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
          #print(paste0("Row strand most similar",row.names(track_row_strands)[maxes],
                       #"highest row overlap",max(get_highest_row)))
        }
        
      }
      #track_row_strands[compare_module,compare_run] <- max(get_highest_row)
    }
    if(max(get_highest_row) < 0.5){
      track_row_strands <- rbind(track_row_strands,rep(0,length(unique_collection)))
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
      track_col_strands <- rbind(track_col_strands,rep(0,length(unique_collection)))
      track_col_strands[nrow(track_col_strands),compare_run] <- 1
      row.names(track_col_strands)[nrow(track_col_strands)] <- paste(compare_run,compare_module,sep=":")
      # add a new module to the comparison set
      col_strands[[length(col_strands)+1]] = compare_col
    }
  }
}

