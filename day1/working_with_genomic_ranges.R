library("rtracklayer")
library("GenomicRanges")

## importing data
## data "100_Morgan_RBiocForAll/CpGislands.Hsapiens.hg38.UCSC.bed" is incluede in the repository
fname <- file.choose()
file.exists(fname)

cpg <- import(fname)
cpg

## work only with the "standard" chromosome 1 - 22 autosomal, X, and Y chromosomes
cpg <- keepStandardChromosomes(cpg, pruning.mode = "coarse")
cpg

head(seqnames(cpg))
head(start(cpg))
head(end(cpg))
head(strand(cpg))

head(cpg$name)

# Extract a vector of width of each CpG island, transform using log10(), draw a histogram
hist(log10(width(cpg)))

# select the CpG islands on chromosomes 1 and 2
subset(cpg,seqnames %in% c("chr1","chr2"))

## Genomic annotations
library("TxDb.Hsapiens.UCSC.hg38.knownGene")

# Extract the coordinates of all transcripts
tx <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
tx

# keep only the standard chromosomes to work with cpg object
tx <- keepStandardChromosomes(tx, pruning.mode = "coarse")
tx

# count the number of CpG islands that overlap each transcript
olaps <- countOverlaps(tx,cpg)
length(olaps)
table(olaps)

# add countOverlaps() to Granges object
tx$cpgOverlaps <- countOverlaps(tx,cpg)
tx

