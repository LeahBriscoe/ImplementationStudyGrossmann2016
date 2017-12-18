# PIPELINE E: GSVA analysis
# Unfinished code

source("https://bioconductor.org/biocLite.R")
#biocLite("reactome.db")
require("reactome.db")
require("GSEABase")
require("GO.db")
require(GSVA)
biocLite("hgu95av2.db")
require("hgu95av2.db")
#biocLite("GO.db")

dataset <- "CleanData1"
export_path <- "~/Documents/MII_Rotation/GrossmannReplication/UnofficialData/"
#export_path <- "/u/home/b/briscoel/scratch/GrossmanLB_Pipeline/"
#idat <- read.csv(paste0(export_path,dataset,"Imaging.csv"),row.names=1,header=TRUE,check.names=FALSE)
edat <- read.csv(paste0(export_path,dataset,"Expression.csv"),row.names=1,header=TRUE,check.names=FALSE)
edat <- as.matrix(edat)

my_eset <- ExpressionSet(assayData=edat)
my_gene_sets <- GeneSet(row.names(my_eset), geneIdType=EntrezIdentifier())
my_gene_set_collection <- GeneSetCollection(my_gene_sets)

dim(edat)
n= 262
p= 21766
nGS <- 100
min.sz <- 10  ## minimum gene set size
max.sz <- 100 ## maximum gene set size
gs <- as.list(sample(min.sz:max.sz, size=nGS, replace=TRUE))
gs <- lapply(gs, function(n, p) sample(1:p, size=n, replace=FALSE), p)

es.max <- gsva(edat, gs, annotation=reactome.db,mx.diff=FALSE, verbose=FALSE, parallel.sz=1)

es.dif <- gsva(edat, gs, mx.diff=TRUE, verbose=FALSE, parallel.sz=1)


my_gsva <- gsva(edat,gset.idx.list=my_gene_set_collection,annotation=reactomeEXTID2PATHID,method='gsva')


gsc <- GeneSetCollection(sample.ExpressionSet[200:250],
                         setType = GOCollection())

my_gsva <- gsva(edat,gsc, annotation=GO.db,mx.diff=FALSE, verbose=FALSE, parallel.sz=1)



reactomePATHID2EXTID(row.names(edat)[1:50])

install.packages("GO")
my_gene_set_collection <- GeneSetCollection(my_gene_sets,idType=EntrezIdentifier())
lst <- head(as.list(reactomeEXTID2PATHID))


mapply
fl <- system.file("extdata", "Broad.xml", package="GSEABase")
gs2 <- getBroadSets(fl)[[1]] # actually, a list of two gene sets

## GeneSet from list of gene identifiers
geneIds <- geneIds(gs2) # any character vector would do
gs3 <- GeneSet(geneIds)
my_gene_sets <- GeneSet(type = "GeneIdentifierType",

