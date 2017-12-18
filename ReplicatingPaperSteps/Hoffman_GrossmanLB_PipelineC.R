# If running shell script "Hoffman_GrossmanLB_Pipeline.sh", pipeline C comes after running the GSEA binary on commandline 
# Must run the bash code after "#Gene Set Enrichment Analysis" before Pipeline C
# PIPELINE-C: Searches for GSEA reports for both negative and positive enrichment for all imaging features and merges
#             all reports into a single NES matrix
# OUTPUT: Single NES Matrix in csv format


# Take in parameters from shell script: 1 parameter: "Data1" or "Data2"
args <- commandArgs(trailingOnly = TRUE)
dataset = args[1]
#dataset= "Data1"
#dir.create(Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
#install.packages(c("biomaRt","annotate","org.Hs.eg.db"), Sys.getenv("R_LIBS_USER"), repos = "http://bioconductor.org/biocLite.R" )
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("biomaRt","annotate","org.Hs.eg.db"))

#gsea_folder <- paste0("/u/home/b/briscoel/scratch/GrossmanLB_Pipeline/",dataset,"GSEA_Analysis/")
gsea_folder <- paste0("~/Documents/MII_Rotation/ModuleRobustness/UnofficialData/",dataset,"EntrezGeneWeightsPerRadiomicFeature/")
#export_path <- "/u/home/b/briscoel/scratch/GrossmanLB_Pipeline/"
export_path <- "~/Documents/MII_Rotation/ModuleRobustness/"


radiomic_folder_names <- dir(path = gsea_folder, pattern="imaging.radiomics.*",full.names = FALSE, recursive=FALSE)
#to be used for colnames
clean_radiomic_names <- gsub("imaging.radiomics.","",radiomic_folder_names)

nes_neg_file_names <- dir(path = gsea_folder,pattern="gsea_report_for_na_neg.*xls",full.names = TRUE, recursive = TRUE)
nes_pos_file_names <- dir(path = gsea_folder,pattern="gsea_report_for_na_pos.*xls",full.names = TRUE, recursive = TRUE)
#print(nes_neg_file_names)

nes_matrix <- matrix(NA,511,1)
for(n_nes in 1:length(nes_neg_file_names)){
  if(n_nes == 1){
    neg_nes <- read.table(nes_neg_file_names[n_nes],sep="\t",header = TRUE,row.names=1)
    pos_nes <- read.table(nes_pos_file_names[n_nes],sep="\t",header = TRUE,row.names=1)
    row.names(nes_matrix) <- c(row.names(neg_nes),row.names(pos_nes))
    nes_matrix[,n_nes] <- c(neg_nes$NES,pos_nes$NES)
  }
  else{
    neg_nes <- read.table(nes_neg_file_names[n_nes],sep="\t",header = TRUE,row.names=1)
    pos_nes <- read.table(nes_pos_file_names[n_nes],sep="\t",header = TRUE,row.names=1)
    temp_nes_list <- c(neg_nes$NES,pos_nes$NES)
    names(temp_nes_list) <- c(row.names(neg_nes),row.names(pos_nes))
    nes_matrix <- cbind(nes_matrix,temp_nes_list[row.names(nes_matrix)])
  }
}

colnames(nes_matrix) <- clean_radiomic_names
write.table(nes_matrix,paste0(export_path,dataset,"_NES_matrix.txt"),sep="\t",quote=FALSE)

