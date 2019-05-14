## RNA-seq analysis
## analyse RNA-seq datafrom the mouse mammary gland, use edgeR package to import, organise, filter and 
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

glMDSPlot(lcpm, labels=paste(group, lane, sep="_"), 
          groups=x$samples[,c(2,5)], launch=FALSE)


title(main="B. Sequencing lanes")
