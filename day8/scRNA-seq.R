## Analysis of single-cell RNA-seq data: dimensionality reduction, clustering, and lineage inference

# The data we use is from a  a scRNA-seq study of stem cell differentiation in the mouse olfactory epithelium 
# (OE) (Fletcher et al. 2017).

suppressPackageStartupMessages({
  # Bioconductor
  library(BiocParallel)
  library(SingleCellExperiment)
  library(clusterExperiment)
  library(scone)
  library(zinbwave)
  library(slingshot)
  # CRAN
  library(gam)
  library(RColorBrewer)
})

set.seed(20)

register(SerialParam())

# counts for all genes in each cell are available as part of the GitHub R package 'fletcher2017data', but it is not 
# available for R 3.6.0, here we just download it from the repository and load them
load('fletcher.rda')
fletcher
colData(fletcher)

# pre-processing
# QC-metric-based sample filtering
data("housekeeping")
hk = rownames(fletcher)[toupper(rownames(fletcher)) %in% housekeeping$V1]
mfilt <- metric_sample_filter(counts(fletcher),
                              nreads = colData(fletcher)$NREADS,
                              ralign = colData(fletcher)$RALIGN,
                              pos_controls = rownames(fletcher) %in% hk,
                              zcut = 3, mixture = FALSE,
                              plot = TRUE)

# simplify to a single logical
mfilt <- !apply(simplify2array(mfilt[!is.na(mfilt)]), 1, any)
filtered <- fletcher[, mfilt]
dim(filtered)
filtered <- makeFilterStats(filtered, filterStats="var", transFun = log1p)
filtered <- filterData(filtered, percentile=1000, filterStats="var")
filtered

publishedClusters <- colData(filtered)[, "publishedClusters"]
col_clus <- c("transparent", "#1B9E77", "antiquewhite2", "cyan", "#E7298A", 
              "#A6CEE3", "#666666", "#E6AB02", "#FFED6F", "darkorchid2", 
              "#B3DE69", "#FF7F00", "#A6761D", "#1F78B4")
names(col_clus) <- sort(unique(publishedClusters))
table(publishedClusters)

# Normalization and dimensionality reduction: ZINB-WaVE
clustered <- zinbwave(filtered, K = 50, X = "~ Batch", residuals = TRUE, normalizedValues = TRUE)

# We already have clustered data in fletcher2017 data package, so we can load it to avoid waiting for the computations
load('clustered.rda')
clustered

# Normalization
assayNames(clustered)
norm <- assay(clustered, "normalizedValues")
norm[1:3,1:3]

# Dimensionality reduction
reducedDimNames(clustered)
W <- reducedDim(clustered, "zinbwave")
W[1:3,1:3]

W <- reducedDim(clustered)
d <- dist(W)
fit <- cmdscale(d, eig = TRUE, k = 2)
plot(fit$points, col = col_clus[as.character(publishedClusters)], main = "",
     pch = 20, xlab = "Component 1", ylab = "Component 2")
legend(x = "topleft", legend = unique(names(col_clus)), cex = .5, fill = unique(col_clus), title = "Sample")

# cell clustering: RSEC
clustered <- RSEC(clustered, k0s = 4:15, alphas = c(0.1),
                  betas = 0.8, reduceMethod="zinbwave",
                  clusterFunction = "hierarchical01", minSizes=1,
                  ncores = NCORES, isCount=FALSE,
                  dendroReduce="zinbwave",
                  subsampleArgs = list(resamp.num=100,
                                       clusterFunction="kmeans",
                                       clusterArgs=list(nstart=10)),
                  verbose=TRUE,
                  consensusProportion = 0.7,
                  mergeMethod = "none", random.seed=424242,
                  consensusMinSize = 10)
clustered
is(clustered, "SingleCellExperiment")
slotNames(clustered)

plotClusters(clustered)
plotCoClustering(clustered)

table(primaryClusterNamed(clustered))

plotBarplot(clustered, legend = FALSE)


clustered <- addClusterings(clustered, colData(clustered)$publishedClusters, 
                            clusterLabel = "publishedClusters")

clusterLegend(clustered)$publishedClusters[, "color"] <- 
  col_clus[clusterLegend(clustered)$publishedClusters[, "name"]]

plotBarplot(clustered, whichClusters=c("makeConsensus", "publishedClusters"),
            xlab = "", legend = FALSE,missingColor="white")

plotClustersTable(clustered, whichClusters=c("makeConsensus","publishedClusters"))

experimentColors <- bigPalette[1:nlevels(colData(clustered)$Experiment)]
batchColors <- bigPalette[1:nlevels(colData(clustered)$Batch)]
metaColors <- list("Experiment" = experimentColors,
                   "Batch" = batchColors)

plotHeatmap(clustered, 
            whichClusters = c("makeConsensus","publishedClusters"), clusterFeaturesData = "all",
            clusterSamplesData = "dendrogramValue", breaks = 0.99,
            colData = c("Batch", "Experiment"),
            clusterLegend = metaColors, annLegend = FALSE, main = "")

plotReducedDims(clustered,whichCluster="primary",reducedDim="zinbwave",pch=20,
                xlab = "Component1", ylab = "Component2",legendTitle="Sample",main="",
                plotUnassigned=FALSE
)

# cell lineage and pseudotime inference: Slingshot
table(data.frame(original = publishedClusters, ours = primaryClusterNamed(clustered)))

pseudoCe <- clustered[,!primaryClusterNamed(clustered) %in% c("-1")]
X <- reducedDim(pseudoCe,type="zinbwave")
mds <- cmdscale(dist(X), eig = TRUE, k = 4)
lineages <- slingshot(mds$points, clusterLabels = primaryClusterNamed(pseudoCe), start.clus = "c1")

colorCl<-convertClusterLegend(pseudoCe,whichCluster="primary",output="matrixColors")[,1]
pairs(lineages, type="lineages", col = colorCl)

pairs(lineages, type="curves", col = colorCl)

lineages

reducedDim(pseudoCe, "MDS") <- mds$points
pseudoCe <- slingshot(pseudoCe, reducedDim = "MDS", start.clus = "c1")
pseudoCe
colData(pseudoCe)

pseudoCeSup <- slingshot(pseudoCe, reducedDim = "MDS", start.clus = "c1",
                         end.clus = c("c3", "c6", "c2"))

# differential expression analysis along lineages
t <- colData(pseudoCe)$slingPseudotime_1
y <- transformData(pseudoCe)
gam.pval <- apply(y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})

topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]

pseudoCe1 <- pseudoCe[,!is.na(t)]
orderSamples(pseudoCe1)<-order(t[!is.na(t)])

plotHeatmap(pseudoCe1[topgenes,], clusterSamplesData = "orderSamplesValue", breaks = .99)

