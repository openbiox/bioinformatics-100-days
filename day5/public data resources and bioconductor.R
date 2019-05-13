required_pkgs = c(
  "TCGAbiolinks", 
  "GEOquery", 
  "GenomicDataCommons",
  "limma",
  "curatedTCGAData",
  "recount",
  "curatedMetagenomicData",
  "phyloseq",
  "HMP16SData",
  "caTools",
  "piano",
  "isa",
  "VennDiagram",
  "downloader",
  "gdata",
  "AnnotationDbi",
  "hgu133a.db",
  "PharmacoGx")
BiocManager::install(required_pkgs)

# GEOquery
library(GEOquery)
gse = getGEO("GSE103512")[[1]]
# convert from ExpressionSet to SummarizedExperiment
library(SummarizedExperiment)
se = as(gse,'SummarizedExperiment')

# examine variables of interest
with(colData(se),table(`cancer.type.ch1`,`normal.ch1`))


# filter gene expression by variance to find most informative genes
sds = apply(assay(se,'exprs'),1,sd)
dat = assay(se,'exprs')[order(sds,decreasing = TRUE)[1:500],]


# perform multidimensional scaling
# make a data frame 
mdsvals = cmdscale(dist(t(dat)))
mdsvals = as.data.frame(mdsvals)
mdsvals$Type=factor(colData(se)[,'cancer.type.ch1'])
mdsvals$Normal = factor(colData(se)[,'normal.ch1'])
head(mdsvals)

library(ggplot2)
ggplot(mdsvals, aes(x=V1,y=V2,shape=Normal,color=Type)) + 
  geom_point( alpha=0.6) + theme(text=element_text(size = 18))



# GenomicDataCommons
library(GenomicDataCommons)

# check network connectivity
GenomicDataCommons::status()

#find data
# build a manifest that can be used to guide the download of raw data
ge_manifest = files() %>%
  filter( ~ cases.project.project_id == 'TCGA-OV' &
            type == 'gene_expression' &
            analysis.workflow_type == 'HTSeq - Counts') %>%
  manifest()

# download data
fnames = lapply(ge_manifest$id[1:20],gdcdata)

# metadata queries
expands = c("diagnoses","annotations",
            "demographic","exposures")
projResults = projects() %>%
  results(size=10)
str(projResults,list.len=5)

names(projResults)

# Querying metadata 
# creating a query, specific subclasses of GDQuery: projects(), cases(), files(), annotations()
pquery = projects()
str(pquery)

# retrieving results
# count() can get records available that satisfy the filter criteria
pcount = count(pquery)
pcount = pquery %>% count()
pcount

# result() fetch actual results
presults = pquery %>% results()

str(presults) # str() is useful for taking a quick glimpse of the data

length(ids(presults))

# results_all() fetch all the available results given a query
presults = pquery %>% results_all()
length(ids(presults))

# include all records
length(ids(presults)) == count(pquery)

# fields and values
default_fields('files')
length(available_fields('files'))
head(available_fields('files'))

qcases = cases()
qcases$fields

qcases = cases() %>% GenomicDataCommons::select(available_fields('cases'))
head(qcases$fields)


# facets and aggregation
# The GDC API offers a feature known as aggregation or faceting
res = files() %>% facet(c('type','data_type')) %>% aggregations()
res$type

# filtering
qfiles = files()
qfiles %>% count()

# filter file results to only "gene_expression" files
qfiles = files() %>% filter(~ type == 'gene_expression')
str(get_filter(qfiles))

# create a filter based on the project
grep('pro',available_fields('files'),value=TRUE)

files() %>% facet('cases.project.project_id') %>% aggregations()

qfiles = files() %>%
  filter( ~ cases.project.project_id == 'TCGA-OV' & type == 'gene_expression')
str(get_filter(qfiles))
qfiles %>% count()

# generate a manifest for buld downloads 
manifest_df = qfiles %>% manifest()
head(manifest_df)

qfiles = files() %>% filter( ~ cases.project.project_id == 'TCGA-OV' &
                               type == 'gene_expression' &
                               analysis.workflow_type == 'HTSeq - Counts')
manifest_df = qfiles %>% manifest()
nrow(manifest_df)


# Datafile access and download
fnames = gdcdata(manifest_df$id[1:2],progress = FALSE)
fnames = gdcdata(manifest_df$id[3:10], access_method = 'client')

# Sequence Read Archive
library(SRAdbV2)

# create a new instance of Omicidx class
oidx = Omicidx$new()

# Queries
query=paste(
  paste0('sample_taxon_id:', 10116),
  'AND experiment_library_strategy:"rna seq"',
  'AND experiment_library_source:transcriptomic',
  'AND experiment_platform:illumina')
z = oidx$search(q=query,entity='full',size=100L)

# fetching results
s = z$scroll()
s

# collating entire result sets
res = s$collate(limit = 1000)
head(res)

# To resus a Scroller, we need reset.
s$reset()
s

# yielding chunks
j = 0
while(s$fetched < 500) {
  res = s$yield()
  # do something interesting with `res` here if you like
  j = j + 1
  message(sprintf('total of %d fetched records, loop iteration # %d', s$fetched, j))
}


## Accessing The Cancer Genome Atlas (TCGA)
#  TCGAbiolinks
library(TCGAbiolinks)
library(SummarizedExperiment)
query <- GDCquery(project = "TCGA-ACC",
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  platform = "Illumina HiSeq", 
                  file.type  = "normalized_results",
                  experimental.strategy = "RNA-Seq",
                  legacy = TRUE)

gdcdir <- file.path("Waldron_PublicData", "GDCdata")
GDCdownload(query, method = "api", files.per.chunk = 10,
            directory = gdcdir)
ACCse <- GDCprepare(query, directory = gdcdir)
ACCse

# curatedTCGAData: Curated Data From The Cancer Genome Atlas as MultiAssayExperiment Objects
library(curatedTCGAData)
library(MultiAssayExperiment)

# By default, the curatedTCGAData() function will only show available datasets, and not download anything.
curatedTCGAData(diseaseCode = "*", assays = "*")
curatedTCGAData(diseaseCode = "ACC")
ACCmae <- curatedTCGAData("ACC", c("RPPAArray", "RNASeq2GeneNorm"), 
                          dry.run=FALSE)
ACCmae

dim(colData(ACCmae))
head(colnames(colData(ACCmae)))

head(metadata(colData(ACCmae))[["subtypes"]])


##  recount: Reproducible RNA-seq Analysis Using recount2
library(recount)
project_info <- abstract_search('GSE32465')

download_study(project_info$project)
load(file.path(project_info$project, 'rse_gene.Rdata'))


## curated*Data packages for standardized cancer transcriptomes
# microbiome data
# curatedMetagenomicData: Curated and processed metagenomic data through ExperimentHub
library(curatedMetagenomicData)
View(data.frame(combined_metadata))

oral <- c("BritoIL_2016.metaphlan_bugs_list.oralcavity",
          "Castro-NallarE_2015.metaphlan_bugs_list.oralcavity")
esl <- curatedMetagenomicData(oral, dryrun = FALSE)

ExpressionSet2phyloseq( esl[[1]], phylogenetictree = TRUE)


#  HMP16SData: 16S rRNA Sequencing Data from the Human Microbiome Project
suppressPackageStartupMessages(library(HMP16SData))
V13()


## Pharmacogenomics
library(PharmacoGx)

# get a list of all the available PharmacoSets in PharmacoGx
availablePSets(saveDir=file.path(".", "Waldron_PublicData"))

# drug sensitivity datasets
psets <- availablePSets(saveDir=file.path(".", "Waldron_PublicData"))
psets[psets[ , "Dataset.Type"] == "sensitivity", ]

# drug perturbation datasets
psets <- availablePSets(saveDir=file.path(".", "Waldron_PublicData"))
psets[psets[ , "Dataset.Type"] == "perturbation", ]

