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
 
 # Differential expression analysis
 # standard DE steps
 dds <- DESeq(dds)
 res <- results(dds) # res contains the results for each gene
 head(res[order(res$pvalue),])
 plotCounts(dds, which.min(res$pvalue),"dex") 
 plotMA(res,ylim = c(-5,5))
 
  # shirinking LFC
 library("apeglm")
 resultsNames(dds)
 res2 <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm") 
 
 par(mfrow=c(1,2))
 plotMA(res, ylim=c(-3,3), main="No shrinkage")
 plotMA(res2, ylim=c(-3,3), main="apeglm")
 
 # minimum effect size
 # perfor a threshold test
 res.lfc <- results(dds, lfcThreshold = 1)
 res.lfc2 <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm", lfcThreshold = 1) 
 par(mfrow=c(1,2))
 plotMA(res.lfc, ylim=c(-5,5), main="No shrinkage, LFC test")
 plotMA(res.lfc2, ylim=c(-5,5), main="apeglm, LFC test", alpha=0.01)
 
 # AnnotationHub
 # Querying AnnotationHub
 library("AnnotationHub")
 ah <- AnnotationHub() 
 display(ah)
 query(ah, c("OrgDb","Homo sapiens"))
 hs <- ah[["AH61777"]]
 hs
 
 # Mapping IDs
 columns(hs)
 table(rownames(res) %in% keys(hs, "ENSEMBL"))
 res$symbol <- mapIds(hs, rownames(res), column="SYMBOL", keytype="ENSEMBL")
 head(res)
 
 # building reports
 # 1. ReportingTools
 library("ReportingTools")
 tmp <- tempdir() # you would instead use a meaningful path here
 rep <- HTMLReport(shortName="airway", title="Airway DGE",
                   basePath=tmp, reportDirectory="report")
 publish(res, rep, dds, n=20, make.plots=TRUE, factor=dds$dex)
 finish(rep)
 
 browseURL(file.path(tmp,"report","airway.html"))
 
 # 2. Glimma
 library("Glimma")
 status <- as.numeric(res$padj < .1)
 anno <- data.frame(GeneID=rownames(res), symbol=res$symbol)
 glMDPlot(res2, status=status, counts=counts(dds,normalized=TRUE),
          groups=dds$dex, transform=FALSE,
          samples=colnames(dds), anno=anno,
          path=tmp, folder="glimma", launch=FALSE)
 browseURL(file.path(tmp,"glimma","MD-Plot.html"))
 
 # integration with ZINB-WaVE
 # simulate single-cell count with splatter
 library("splatter")
 params <- newSplatParams()
 params <- setParam(params, "de.facLoc", 1) 
 params <- setParam(params, "de.facScale", .25)
 params <- setParam(params, "dropout.type", "experiment")
 params <- setParam(params, "dropout.mid", 3)
 
 set.seed(1)
 sim <- splatSimulate(params, group.prob=c(.5,.5), method="groups")
 
 plot(log10(rowMeans(assays(sim)[["TrueCounts"]])),
      rowMeans(assays(sim)[["Dropout"]]))
 
 rowData(sim)$log2FC <- with(rowData(sim), log2(DEFacGroup2/DEFacGroup1))
 
 rowData(sim)$trueDisp <- rowMeans(assays(sim)[["BCV"]])^2
 gridlines <- c(1e-2,1e-1,1); cols <- c("blue","red","darkgreen")
 with(rowData(sim)[rowData(sim)$GeneMean> 1,],
      plot(GeneMean, trueDisp, log="xy", xlim=c(1,300), ylim=c(.01,5)))
 abline(h=gridlines, col=cols)
 text(300, gridlines, labels=gridlines, col=cols, pos=3)
 
 # model zeros with zinbwave 
 library(zinbwave)
 keep <- rowSums(counts(sim) >= 5) >= 25
 table(keep)
 
 zinb <- sim[keep,]
 zinb$condition <- factor(zinb$Group)
 
 nms <- c("counts", setdiff(assayNames(zinb), "counts"))
 assays(zinb) <- assays(zinb)[nms]
 
 zinb <- zinbwave(zinb, K=0, BPPARAM=SerialParam(), epsilon=1e12)
 
 # model non-zeros with DESeq2
 zdds <- DESeqDataSet(zinb, design=~condition)
 zdds <- DESeq(zdds, test="LRT", reduced=~1,
               sfType="poscounts", minmu=1e-6, minRep=Inf)
 
 # plot dispersion estimates
 plotDispEsts(zdds)
 
 keepForDispTrend <- rowSums(counts(zdds) >= 10) >= 25
 zdds2 <- estimateDispersionsFit(zdds[keepForDispTrend,])
 plotDispEsts(zdds2)
 
 dispersionFunction(zdds) <- dispersionFunction(zdds2)
 zdds <- estimateDispersionsMAP(zdds)
 
 zdds <- nbinomLRT(zdds, reduced=~1, minmu=1e-6)
 
 # evaluation against truth 
 with(mcols(zdds), plot(trueDisp, dispMAP, log="xy"))
 abline(0,1,col="red")
 zres <- results(zdds, independentFiltering=FALSE)
 plot(mcols(zdds)$log2FC, zres$log2FoldChange, ylim=c(-4,4)); abline(0,1,col="red")
 
 ncts <- counts(zdds, normalized=TRUE)
 simple.lfc <- log2(rowMeans(ncts[,zdds$condition == "Group2"])/
                      rowMeans(ncts[,zdds$condition == "Group1"]))
 plot(mcols(zdds)$log2FC, simple.lfc, ylim=c(-4,4)); abline(0,1,col="red")
 
 tab <- table(sig=zres$padj < .05, DE.status=mcols(zdds)$log2FC != 0)
 tab
 
 round(prop.table(tab, 1), 3)
 
 