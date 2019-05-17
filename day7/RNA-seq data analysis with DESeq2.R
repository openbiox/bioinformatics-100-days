## RNA-seq data analysis with DESeq2

 library(tximportData)
 
 # import Salmon transcript abundances squantification files using the data in the tximportData package
 dir <- system.file("extdata", package="tximportData") # function system.file can be used to find out where on your 
                                                       # computer the files from a package have been installed
 list.files(dir)
 list.files(file.path(dir,"salmon"))
 
 samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
 samples
 files <- file.path(dir, "salmon", samples$run, "quant.sf.gz")
 names(files) <- paste0("sample",1:6)
 all(file.exists(files))

 ## Mapping transcripts to genes
 library("TxDb.Hsapiens.UCSC.hg38.knownGene")
 txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
 k <- keys(txdb, keytype="TXNAME")
 tx2gene <- select(txdb, k, "GENEID", "TXNAME")
 
 library("readr")
 tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
 head(tx2gene)

 ## tximport command
 library("tximport")
 library("jsonlite")
 library("readr")
 txi <- tximport(files, type="salmon", tx2gene=tx2gene) # txi object is a list of matrices
 
 names(txi)
 txi$counts[1:3,1:3]
 txi$length[1:3,1:3]
 txi$abundance[1:3,1:3]
 txi$countsFromAbundance
 
 library("DESeq2")
 dds <- DESeqDataSetFromTximport(txi, samples, ~1) # design formula of ~1 means we can only fit an intercept term
 dds$center
 dds$pop
 
 ## Exploratory data analysis
 # simple EDA
 library("airway")
 data("airway") 
 airway$dex <- relevel(airway$dex, "untrt") 
 airway$dex

 round( colSums(assay(airway)) / 1e6, 1 )
 
 colData(airway)
 table(airway$cell)
 table(airway$dex)
 
 dds <- DESeqDataSet(airway, design = ~ cell + dex)
 keep <- rowSums(counts(dds) >= 5) >= 4 # minimal filtering
 table(keep)
 dds <- dds[keep,]
 boxplot(log10(counts(dds)+1))
 
 dds <- estimateSizeFactors(dds)
 boxplot(log10(counts(dds,normalized=TRUE)+1))
 
 ## Data transformation for EDA
 vsd <- vst(dds)
 class(vsd) 

 assay(vsd)[1:3,1:3]
 
 # pricipal components plot
 plotPCA(vsd, "dex")
 library("ggplot2")
 pcaData <- plotPCA(vsd, intgroup = c( "dex", "cell"), returnData = TRUE)
 percentVar <- round(100 * attr(pcaData, "percentVar"))
 ggplot(pcaData, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
   geom_point(size =3) +
   xlab(paste0("PC1: ", percentVar[1], "% variance")) +
   ylab(paste0("PC2: ", percentVar[2], "% variance")) +
   coord_fixed()