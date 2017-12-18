# ImplementationStudyGrossmann2016
Replicating computation of modules of molecular pathways and radiomic (imaging) features from lung CT radiomics data published in a 2016 eLife publication by Grossmann et al. "Defining the biological basis of radiomic phenotypes in lung cancer" (1). Assessing stability of the Iterative Signature Algorithm for producing modules. 

# Description of Folders
## Replicating Paper Steps
Using the publication as a guideline, I implemented the steps used to create modules, trying to copy these methods as closely as possible. I used the same version of C2 Reactome Gene Sets v4.0, which can be found at http://software.broadinstitute.org/gsea/downloads_archive.jsp. You will need to download the version 4.0 zip file to access the Reactome Set. Also, you will need to download the jar file for javaGSEA at http://software.broadinstitute.org/gsea/downloads.jsp.

## StudyingISA
Scripts to understand ISA output and compare parameters.

Header file in each script describes it's purpose

## ExampleModuleTxtFiles
Example of what is output by ReplicatingPaperSteps>Hoffman_GrossmannLB_PipelineD.R. The txt files contain lists of imaging features and pathways in one module.

## Unofficial Data
The xlsx files contain the Imaging, Expression and Clinical data provided by the publication. I have created my own csv files for each of the Imaging, Expression and Clinical data. Other files are my own Normalized Enrichment Score (NES) matrices. Some files were too large to upload to GitHub. Please check https://drive.google.com/open?id=1RS-jIIub4yjYYWgJ4J9AGU5jouRyps7Y for the full set.

# Description of Miscellaneous Files
## GrossmannsModules.csv
Using the PDF of module heatmaps provided in the supplement to the publication, I extracted text to create a reference file listing out the imaging features and pathways in each of the 13 modules described in the publication.
## CompareModules.py
Using the .txt files describing the modules I produce (from ReplicatingPaperSteps>Hoffman_GrossmannLB_PipelineD.R), I compare to the modules in the publication (in GrossmannsModules.csv). Running python CompareModules.py will print out the highest similarity (highest overlap) between each of my modules with Grossmann's modules - also printing out which module from 1 to 13 was most similar to my modules.



1.Grossmann, P. et al. Defining the biological basis of radiomic phenotypes in lung cancer. eLife 6, (2017).
