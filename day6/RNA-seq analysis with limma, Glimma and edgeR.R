## RNA-seq analysis
## analyse RNA-seq data from the mouse mammary gland, use edgeR package to import, organise, filter and 
## normalize the data, yse limma package with its voom method, linear modelling and empirical Bayes moderation
## to asses differential expression and graphical representations
## Glimma package enables interactive exploration of the results so that individual samples and genes can be 
## examined by the user

# install RNAseq123 workflow
source("https://bioconductor.org/biocLite.R")
biocLite("RNAseq123")

suppressPackageStartupMessages({
  library(limma)
  library(Glimma)
  library(edgeR)
  library(Mus.musculus)
})

# Data packaging
dir.create("Law_RNAseq123")
setwd("Law_RNAseq123")

url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"  
utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb")   
utils::untar("GSE63310_RAW.tar", exdir = ".")  
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
           "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
           "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")
for(i in paste(files, ".gz", sep=""))  
  R.utils::gunzip(i, overwrite=TRUE) 

# reading in count-data
read.delim(file.path(files[1]), nrow=5)

# edgeR offers a convenient way to read the nine text files separately and combine them into a matrix of 
## counts in one step
x <- readDGE(file.path(files),columns = c(1,3))
class(x)
dim(x)
# if the counts from all samples were stored in a single file, the data can be read into R and then converted into a 
# DGEList-object using the DGEList function

# Organising sample information
# sample-level information related to the experimwntal design needs to be associated with the columns of the counts matrix
# DGEList-object contains a samples data frame that stores both cell type and batch information
samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
samplenames
colnames(x) <- samplenames
group <- as.factor(c("LP","ML","Basal","Basal","ML","LP","Basal","ML","LP"))
x$samples$group <- group
lane <- as.factor((rep(c("L004","L006","L008"),c(3,4,2))))
x$samples$lane <- lane
x$samples

# Organising gene annotations
# second data frame named genes in the DGEList-object is used to store gene-level information associated with rows of 
# the counts matrix
geneid <- rownames(x)
genes <- select(Mus.musculus, keys = geneid, columns = c("SYMBOL", "TXCHROM"),keytype = 'ENTREZID')
head(genes)

# keep only the first occurence of each gene ID
genes <- genes[!duplicated(genes$ENTREZID),]

x$genes <- genes
x

## Data pre-processing
# transformations from the raw-scale
# it's common practice to transform raw counts onto a scale that accounts for such library size differences
# popular transformations include: 
# counts per millon (CPM),
# log2-counts per millon (log-CPM)
# reads per kilobase of transcript per million (RPKM)
# fragments per kilobase of transcript permillion (FPKM)
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

# removing genes that are lowly expressed
table(rowSums(x$counts==0)==9)

# set the expression cutoff as 1 to separate expressed genes from unexpressed genes
keep.exprs <- rowSums(cpm>1)>=3
x <- x[keep.exprs,,keep.lib.sizes = FALSE]
dim(x)


library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

# normalizing gene expression distributions
# normalization by the method of trimmed mean of M-values (TMM)
x <- calcNormFactors(x,method = "TMM")
x$ samples$norm.factors

# to give a better visual representation of the effects of normalization, the data was duplicated then adjusted so
# that the counts of the first sample are reduced to 5% of their original values, and in the second sample they are
# inflated to be 5-times larger
x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5

# expression distribution of samples for unnormalized and normalized data
par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
x2 <- calcNormFactors(x2)  
x2$samples$norm.factors
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")

# unsupervised clustering of samples
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")

glMDSPlot(lcpm, labels=paste(group, lane, sep="_"), 
          groups=x$samples[,c(2,5)], launch=FALSE)


## Differential expression analysis 
# creating a design matrix and contrasts
# set up a design matrix with cell population and sequencing lane information
design <- model.matrix(~0+group+lane)
# model.matrix creates a design matrix, e.g by expanding factors to a set of dummy variables
colnames(design) <- gsub("group", "", colnames(design))
design

# contrasts for pairwise comparisons between cell populations
contr.matrix <- makeContrasts(
  BasalvsLP = Basal - LP,
  BasalvsML = Basal - ML,
  LPvsML = LP - ML,
  levels = colnames(design)
)
contr.matrix

# removing heteroscedascity from count data 
par(mfrow = c(1,2)) # par can be used to set or query graphical parameters
                    # mfrow: a vector of the form c(nr, nc). Subsequent figures will be drawn in an nr-by-nc array on
                    # the device by rows (mfrow)

v <- voom(x,design,plot = TRUE) # voom converts raw counts to log-CPM values by automatically extracting library sizes 
                                # and normalization factors from x itself
v # voom-plot provides a visual check on the level of filtering performed upstream

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main = "Final model: Mean-Variance trend")

# fitting linear models for comparisons of interest
# examining the number of DE genes
summary(decideTests(efit)) # significance is defined using an adjusted p-value cutoff that is set at 5% by default

# treat method can be used to calculate p-values from empirical Bayes moderated t-statistics with a minumum log-FC requirement
tfit <- treat(vfit, lfc = 1)
dt <- decideTests(tfit) # genes that are differentially expressed in multiple comparisons can be extracted using the results
                        # from decideTests
summary(dt)

de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)
head(tfit$genes$SYMBOL[de.common], n = 20)
vennDiagram(dt[,1:2], circle.col = c("turquoise","salmon"))

write.fit(tfit, dt, file = "results.txt")


## examining individual DE genes from top to bottom 
basal.vs.lp <- topTreat(tfit, coef = 1, n = Inf)# By default topTreat arranges genes from smallest to largest adjusted p-value 
                                                # with associated gene information, log-FC, average log-CPM, moderated t-statistic, 
                                                # raw and adjusted p-value for each gene
basal.vs.ml <- topTreat(tfit, coef = 2, n = Inf)
head(basal.vs.lp)
head(basal.vs.ml)

# useful graphical representations of differential expression results
# mean-difference plots, display log-FCs from the linear model fit against the average log-CPM values
plotMD(tfit, column = 1, status = dt[,1], main = colnames(tfit)[1], xlim = c(-8, 13))

glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         side.main="ENTREZID", counts=x$counts, groups=group, launch=FALSE)

# heatmap for the top 100 DE genes (ranked by adjusted p-value)
library(gplots)
basal.vs.lp.topgenes <- basal.vs.lp$ENTREZID[1:100]
i <- which(v$genes$ENTREZID %in% basal.vs.lp.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(v$E[i,], scale="row",
          labRow=v$genes$SYMBOL[i], labCol=group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")





