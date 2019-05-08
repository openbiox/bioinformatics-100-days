## working with summarized experimental data
library("airway")
data(airway)

# object construction
colData <- as.data.frame(colData(airway))
counts <- as.data.frame(assay(airway))

# creatins a SummarizedExperiment object
library("SummarizedExperiment")
se <- SummarizedExperiment(assay = counts, colData = colData)
se

subset(se, , dex == "trt")

# calculate the library size
colSums(assay(se))

se$lib.size <- colSums(assay(se))
colData(se)

## downstream analysis
library('DESeq2')

dds <- DESeqDataSet(se, design = ~ cell + dex)
dds

dds <- DESeq(dds)

results(dds)