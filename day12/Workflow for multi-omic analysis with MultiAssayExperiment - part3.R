# Workflow for multi-omics analysis with MultiAssayExperiment - part3

# The Cancer Genome Atlas (TCGA) as MultiAssayExperiment objects
library(curatedTCGAData)
curatedTCGAData("ACC")
suppressMessages({
  acc <- curatedTCGAData("ACC",
                         assays = c("miRNASeqGene", "RPPAArray", "Mutation", "RNASeq2GeneNorm", "CNVSNP"),
                         dry.run = FALSE)
})

dim(colData(acc))
tail(colnames(colData(acc)), 10)

# Utilities for TCGA
library(TCGAutils)
# simplifyTCGA function creates a more manageable MultiAssayExperiment object and using RangedSummarizedExperiment assays where possible
(simpa <- TCGAutils::simplifyTCGA(acc))

# sample types in the data
# sampleTables function gives you an overview of samples in each assay
sampleTables(acc)
head(sampleTypes)

# Curated molecular subtypes
getSubtypeMap(acc)
head(colData(acc)$Histology)

# converting TCGA UUIDs to barcodes and back
UUIDtoBarcode()
barcodeToUUID()
UUIDtoUUID()
filenameToBarcode()

# Plotting, correlation, and other analyses
upsetSamples(miniACC)

###############################

# Kaplan-meier plot stratified by pathology_N_stage

# the colData provides clinical data for things like a Kaplan-Meier plot for overall survival stratified by nodal stage
Surv(miniACC$days_to_death, miniACC$vital_status)

# remove any patients missing overall survival information
miniACCsurv <- miniACC[, complete.cases(miniACC$days_to_death, miniACC$vital_status)]

fit <- survfit(Surv(days_to_death, vital_status) ~ pathology_N_stage, data = colData(miniACCsurv))
ggsurvplot(fit, data = colData(miniACCsurv), risk.table = TRUE)

###############################

# Multivariate Cox regression including RNA-seq, copy number, and pathology
# choose the EZH2 gene for demonstration
wideacc = wideFormat(miniACC["EZH2", , ], 
                     colDataCols = c("vital_status", "days_to_death", "pathology_N_stage"))
wideacc$y = Surv(wideacc$days_to_death, wideacc$vital_status)
head(wideacc)

# perform a multivariate Cox regression with EZH2 copy number (gistict), log2-transformed EZH2 expression (RNASeq2GeneNorm), and nodal
# status (pathology_N_stage) as predictors
coxph(Surv(days_to_death, vital_status) ~ gistict_EZH2 + log2(RNASeq2GeneNorm_EZH2) + pathology_N_stage,
      data=wideacc)


# correlation between RNA-seq and copy number

# narrow down miniACC to only the assays needed
subacc <- miniACC[, , c("RNASeq2GeneNorm", "gistict")]
# Align the rows and columns, keeping only samples with both assays available
subacc <- intersectColumns(subacc)
subacc <- intersectRows(subacc)

subacc.list <- assays(subacc)
subacc.list[[1]] <- log2(subacc.list[[1]] + 1 )
subacc.list <- lapply(subacc.list, t)

corres <- cor(subacc.list[[1]], subacc.list[[2]])

hist(diag(corres))

hist(corres[upper.tri(corres)])

## for the gene with highest correlation to copy number, make a box plot of log2 expression against copy number
which.max(diag(corres))
df <- wideFormat(subacc["EIF4E", , ])
head(df)
boxplot(RNASeq2GeneNorm_EIF4E ~ gistict_EIF4E,
        data = df, varwidth = TRUE,
        xlab = "GISTIC Relative Copy Number Call",
        ylab = "RNA-seq counts")

## Identifying correlated principal components
# Performing PCA of each of the five assays, using samples available on each assay, log-transforming RNA-seq data first.
# Using the first 10 components, calculate Pearson correlaiton between all scores and plot these correlations as a heatmap to identify
# correlated components across assays.
getLoadings <- function(x, ncomp=10, dolog=FALSE, center=TRUE, scale.=TRUE){
  if(dolog){
    x <- log2(x + 1)
  }
  pc = prcomp(x, center=center, scale.=scale.)
  return(t(pc$rotation[, 1:10]))
}

miniACC2 <- intersectColumns(miniACC)
miniACC2 <- c(miniACC2, rnaseqPCA=getLoadings(assays(miniACC2)[[1]], dolog=TRUE), mapFrom=1L)
miniACC2 <- c(miniACC2, gistictPCA=getLoadings(assays(miniACC2)[[2]], center=FALSE, scale.=FALSE), mapFrom=2L)
miniACC2 <- c(miniACC2, proteinPCA=getLoadings(assays(miniACC2)[[3]]), mapFrom=3L)
miniACC2 <- c(miniACC2, mutationsPCA=getLoadings(assays(miniACC2)[[4]], center=FALSE, scale.=FALSE), mapFrom=4L)
miniACC2 <- c(miniACC2, miRNAPCA=getLoadings(assays(miniACC2)[[5]]), mapFrom=5L)

miniACC2 <- miniACC2[, , 6:10]
experiments(miniACC2)

df <- wideFormat(miniACC2)[, -1]
mycors <- cor(as.matrix(df))
mycors <- abs(mycors)
diag(mycors) <- NA

has.high.cor <- apply(mycors, 2, max, na.rm=TRUE) > 0.5
mycors <- mycors[has.high.cor, has.high.cor]
pheatmap::pheatmap(mycors)

## Annotating with ranges
getrr <- function(identifiers, EnsDbFilterFunc=AnnotationFilter::SymbolFilter) {
  edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  afl <- AnnotationFilterList(
    EnsDbFilterFunc(identifiers),
    SeqNameFilter(c(1:21, "X", "Y")),
    TxBiotypeFilter("protein_coding"))
  gr <- genes(edb, filter=afl)
  grl <- split(gr, factor(identifiers))
  grl <- grl[match(identifiers, names(grl))]
  stopifnot(identical(names(grl), identifiers))
  return(grl)
}

getrr(rownames(miniACC)[[1]])

rseACC <- miniACC
withRSE <- c(1:3, 5)
for (i in withRSE){
  rowRanges(rseACC[[i]]) <- getrr(rownames(rseACC[[i]]))
}

experiments(rseACC)

# subsetting by ranges
rseACC[GRanges(seqnames="1:1-1e9"), , withRSE]
