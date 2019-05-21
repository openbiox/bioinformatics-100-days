## Functional enrichment analysis of high-throughput omics data

# test whether known biological functions or processes are over-represented (=enriched) in an experimentally-derived gene 
# list 
# Example: Transcriptomic study, in which 12,671 genes have been tested for differential expression between two sample 
# condition and 529 genes were found DE (differentially expressed).
# Among the DE genes, 28 are annotated to a specific functional gene set, which contains in total 170 genes.
deTable <- matrix(c(28, 142, 501, 12000),
                  nrow = 2,
                  dimnames = list(c("DE", "Not.DE"),
                                  c("In.gene.set","Not.in.gene.set")))
deTable

# The overlap of 28 genes can be assessed based on the hypergeometric distribution. This corresponds to a one-sided 
# version of Fisher's exact test, yielding here a significant enrichment.
fisher.test(deTable, alternative = "greater")

# Gene expression-based enrichment analysis
suppressPackageStartupMessages(library(EnrichmentBrowser))

# microarray data
# Here we consider expression measurements of patients with acute lymphoblastic leukemia. A frequent chromosomal defect 
# found among these patients is a translocation, in which parts of chromosome 9 and 22 swap places. This results in the 
# oncogenic fusion gene BCR/ABL created by ABL1 gene on chromosome 9 to a part of the BCR gene on chromosome 22.
library(ALL)
data(ALL)
# select B-cell ALL patients with and without the BCR/ABL fusion
ind.bs <- grep("^B",ALL$BT)
ind.mut <- which(ALL$mol.biol %in% c("BCR/ABL", "NEG"))
sset <- intersect(ind.bs, ind.mut)
all.eset <- ALL[,sset]

dim(all.eset)
exprs(all.eset)[1:4,1:4]

# we often have more than one probe per gene, we compute gene expression values as the average of the corresponding probe values
allSE <- probe2gene(all.eset)
head(names(allSE))

# RNA-seq data
# transcriptome profiles of four primary human airway smooth muscle cell lines in two conditions: control and treatment with 
# dexamethasone
library(airway)
data(airway)

# only keep genes that are annotated to an ENSEMBL gene ID
airSE <- airway[grep("^ENSG", names(airway)), ]
dim(airSE)

assay(airSE)[1:4,1:4]

# Differential expression analysis
# GROUP defines the sample groups being contrasted
# BLOCK defines paired samples or sample blocks

# For the ALL dataset, the GROUP variable indicates whether the BCR-ABL gene fusion is present(1) or not (0)
allSE$GROUP <- ifelse(allSE$mol.biol == "BCR/ABL", 1, 0)
table(allSE$GROUP)

# For the airway data set, it indicates whether the cell lines have been treated with dexamethasone (1) or not (0)
airSE$GROUP <- ifelse(colData(airway)$dex == "trt", 1, 0)
table(airSE$GROUP)

# For the airway dataset, the sample blocks correspond to the four different cell lines
airSE$BLOCK <- airway$cell
table(airSE$BLOCK)

# For microarrya data, the EnrichmentBrowser::deAna function carries out differential expression analysis based on the functionality
# from the limma package. Resulting log2 fold changes and t-test derived p-values for each gene are sppended to the rowData slot.
allSE <- deAna(allSE)
rowData(allSE,use.names = TRUE)

# For RNA-seq data, the deAna function can be used to carry out differential expression analysis between the two groups either based
# on functionality from limma, or alternatively, the frequently used edgeR or DESeq2 package.
airSE <- deAna(airSE, de.method = "edgeR")

rowData(airSE, use.names = TRUE)


## Gene sets
# check whether pre-defined sets of genes that are known to work together
kegg.gs <- getGenesets(org = "hsa", db = "kegg") # download all KEGG pathways for Homo Sapien as gene sets

go.gs <- getGenesets(org = "hsa", db = "go", go.onto = "BP", go.mode = "GO.db") # retrieve GO terms of biological process

# if provided a file, function 'getGenesets' parses user-defined gene sets from GMT file format
# here we use this for reading a list of already downloaded KEGG gene sets
data.dir <- system.file("extdata", package = "EnrichmentBrowser")
gmt.file <- file.path(data.dir, "hsa_kegg_gs.gmt")
hsa.gs <- getGenesets(gmt.file)
hsa.gs[1:2]

## GO/KEGG overrepresentation analysis (ORA)
# ORA is a basic and frequently used method of gene set analysis methods, it tests the overlap between DE genes and
# genes in a gene set based on the hypergeometric distribution
ora.all <- sbea(method = "ora", se = allSE, gs = hsa.gs, perm = 0, alpha = 0.2) # we choose a significanse level alpha = 0.2
gsRanking(ora.all)

eaBrowse(ora.all)

airSE <- idMap(airSE, org = "hsa", from = "ENSEMBL", to = "ENTREZID") # map the airway dataset to Entrez IDs
ora.air <- sbea(method = "ora", se = airSE, gs = hsa.gs, perm = 0)
gsRanking(ora.air)

## Functional class scoring & permutation testing
gsea.all <- sbea(method = "gsea", se = allSE, gs = hsa.gs, perm = 1000)
gsRanking(gsea.all)

# adapted version of GSEA, allows incorporation of limma/voom, edgeR, or DESeq2
gsea.air <- sbea(method = "gsea", se = airSE, gs = hsa.gs, perm = 100)

# use rotation instead of permutation
roast.air <- sbea(method = "roast", se = airSE, gs = hsa.gs)
gsRanking(roast.air)

# additional methods
sbeaMethods()

## Network-based enrichment analysis
# having found gene sets that show enrichment for differential expression, we are now interested whether these findings can be
# supported by known regulatory interactions, e.g: whehter transcription factors and their target genes are expressed in accordance
# to the connecting regulations (activation/inhibition)

# compile a network from regulations in pathway database
hsa.grn <- compileGRN(org = "hsa", db = "kegg")
head(hsa.grn)

# signaling pathway impact analysis (SPIA)
# evaluate whether expression changes are propagated across the pathway topology in combination with ORA
spia.all <- nbea(method = "spia", se = allSE, gs = hsa.gs, grn = hsa.grn, alpha = 0.2)
gsRanking(spia.all)

# gene graph enrichment analysis (GGEA)
ggea.all <- nbea(method = "ggea", se = allSE, gs = hsa.gs, grn = hsa.grn)
gsRanking(ggea.all)

nbeaMethods()
