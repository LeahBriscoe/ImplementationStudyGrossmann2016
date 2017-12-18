#!/bin/bash
# Run R script to get ranks, calculate GSEA NES
. /u/local/Modules/default/init/modules.sh
module load R/3.2.3
echo $1
# Run Pipeline A: Spearman Correlation and Weighting
#R CMD BATCH /u/home/b/briscoel/scratch/GrossmanLB_Pipeline/Hoffman_GrossmanLB_PipelineA.R
# Run Pipeline B: Create Ranked List
# Run Pipeline B_Bicluster: Create Ranked List
#R CMD BATCH /u/home/b/briscoel/scratch/GrossmanLB_Pipeline/Hoffman_GrossmanLB_PipelineB.R

#Gene Set Enrichment Analysis
#module load java/1.8.0_111
# for filename in /u/home/b/briscoel/scratch/GrossmanLB_Pipeline/Data1EntrezGeneWeightsPerRadiomicFeature/*.rnk
# do
#    filenamu=${filename##*/}
#    java -cp /u/home/b/briscoel/scratch/GrossmanLB_Pipeline/gsea-3.0.jar \
#    -Xmx512m xtools.gsea.GseaPreranked \
#    -gmx /u/home/b/briscoel/scratch/GrossmanLB_Pipeline/c2.cp.reactome.v4.0.entrez.gmt \
#    -norm meandiv -nperm 1000 \
#    -rnk $filename \
#    -scoring_scheme weighted -rpt_label my_analysis -create_svgs false -make_sets true -plot_top_x 1 \
#    -rnd_seed 987654321 -set_max 500 -set_min 15 -zip_report false \
#    -out /u/home/b/briscoel/scratch/GrossmanLB_Pipeline/Data1GSEA_Analysis/${filenamu%.*} -gui false;
# done
Rscript /u/home/b/briscoel/scratch/GrossmanLB_Pipeline/Hoffman_GrossmanLB_PipelineC.R $1 > /u/home/b/briscoel/scratch/GrossmanLB_Pipeline/PipelineC.Rout 


# for filename in /u/home/b/briscoel/scratch/GrossmanLB_Pipeline/Data1EntrezGeneWeightsPerRadiomicFeature/*.rnk
# do
#    filenamu=${filename##*/}
#    if [ ! -d /u/home/b/briscoel/scratch/GrossmanLB_Pipeline/Data1GSEA_Analysis/${filenamu%.*} ]; then
# 	   echo ${filenamu%.*}
# 	   java -cp /u/home/b/briscoel/scratch/GrossmanLB_Pipeline/gsea-3.0.jar \
# 	   -Xmx512m xtools.gsea.GseaPreranked \
# 	   -gmx /u/home/b/briscoel/scratch/GrossmanLB_Pipeline/c2.cp.reactome.v4.0.entrez.gmt \
# 	   -norm meandiv -nperm 1000 \
# 	   -rnk $filename \
# 	   -scoring_scheme weighted -rpt_label my_analysis -create_svgs false -make_sets true -plot_top_x 1 \
# 	   -rnd_seed 987654321 -set_max 500 -set_min 15 -zip_report false \
# 	   -out /u/home/b/briscoel/scratch/GrossmanLB_Pipeline/Data1GSEA_Analysis/${filenamu%.*} -gui false;
#    fi
# done