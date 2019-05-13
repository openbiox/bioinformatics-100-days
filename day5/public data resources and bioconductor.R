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
se = as(gse,'SummarisedExperiment')

with(colData(se),table('cancer.type.ch1','normal.ch1'))

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
# creating a query
pquery = projects()
str(query)

# retrieving results
pcount = count(pquery)
pcount = pquery %>% count()
pcount

presults = pquery %>% results()

str(presults)

length(ids(presults))
length(ids(presults))
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
qcases = cases() %>% GenomicDataCommons::select(available_fields('cases'))
head(qcases$fields)

# filtering
qfiles = files()
qfiles %>% count()

qfiles = files() %>% filter(~ type == 'gene_expression')
str(get_filter(qfiles))

grep('pro',available_fields('files'),value=TRUE)

files() %>% facet('cases.project.project_id') %>% aggregations()

qfiles = files() %>%
  filter( ~ cases.project.project_id == 'TCGA-OV' & type == 'gene_expression')
str(get_filter(qfiles))

manifest_df = qfiles %>% manifest()
head(manifest_df)

qfiles = files() %>% filter( ~ cases.project.project_id == 'TCGA-OV' &
                               type == 'gene_expression' &
                               analysis.workflow_type == 'HTSeq - Counts')
manifest_df = qfiles %>% manifest()
nrow(manifest_df)


