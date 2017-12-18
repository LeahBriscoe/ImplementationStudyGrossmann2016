#!/bin/bash
# Run R script to get ranks, calculate GSEA NES
. /u/local/Modules/default/init/modules.sh
module load R/3.2.3

Rscript /u/home/b/briscoel/scratch/GrossmanLB_Pipeline/ImplementationStudy/Hoffman_ISAonGenes.R > /u/home/b/briscoel/scratch/GrossmanLB_Pipeline/ImplementationStudy/ISAonGenes.Rout 
