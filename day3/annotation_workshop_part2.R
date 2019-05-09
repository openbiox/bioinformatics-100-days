# BSgenome packages
# contain sequence information for a given species/build
library(BSgenome)

# get a list of available genomes
head(available.genomes())

# load and inspect a BSgenome package
library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens

# main accessor is "getSeq"
# get data by sequence
getSeq(Hsapiens,"chr1")

# passing in a Granges object, to get just a region
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
gns <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
gns
getSeq(Hsapiens,gns["5467",])

# AnnotationHub
# allow us to query and download many different annotation objects, without having to explicitly install them
library(AnnotationHub)
hub <- AnnotationHub()
hub

# querying AnnotationHub (useful queries are based on Data provider, Data class, Species, Data source)
names(mcols(hub))

# AnnotationHub Data providers
unique(hub$dataprovider)

# AnnotationHub Data classes
unique(hub$rdataclass)

# AnnotationHub species
head(unique(hub$species))

# AnnotationHub Data sources
unique(hub$sourcetype)

# AnnotationHub query
qry <- query(hub,c("granges","homo sapiens","ensembl"))
qry

qry$sourceurl

# selecting AnnotationHub resource
whatIwant <- qry[["AH50377"]]

# convert it to a TxDb format
GRCh38TxDb <- makeTxDbFromGRanges(whatIwant)
GRCh38TxDb

# biomaRt
# biomaRt package allows queries to an Ensembl Biomart server
library(biomaRt)
listMarts()

# check for the available data sets on a particular server
mart <- useMart("ENSEMBL_MART_ENSEMBL")
head(listDatasets(mart))

# biomaRt queries
mart <- useMart("ENSEMBL_MART_ENSEMBL","hsapiens_gene_ensembl") #set up the mart object
atrib <- listAttributes(mart)
filts <- listFilters(mart)
head(atrib)
head(filts)

afyids <- c("1000_at","1001_at","1002_f_at","1007_s_at")
getBM(c("affy_hg_u95av2","hgnc_symbol"),c("affy_hg_u95av2"),afyids,mart)
