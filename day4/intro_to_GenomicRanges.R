# Solving common bioinformatic challenges using GenomicRanges
# 2019-5-10

# Constructing a Granges object from data.frame
suppressPackageStartupMessages({
  library(BiocStyle)
  library(GenomicRanges)
})

df <- data.frame(
  seqnames = rep(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
  start = c(101, 105, 125, 132, 134, 152, 153, 160, 166, 170),
  end = c(104, 120, 133, 132, 155, 154, 159, 166, 171, 190),
  strand = rep(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  score = 1:10,
  GC = seq(1, 0, length=10),
  row.names = head(letters, 10))

gr <- makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)


# Basic manipulation of GRanges objects
# seqnames(), ranges(), strand() accessor funtions extract the components of the genomic coordinates
seqnames(gr)
ranges(gr)
strand(gr)

# granges() function extracts genomic ranges without corresponding metadata
granges(gr)

# the start(), end(), width(), and range functions extract basic interval characteristic
start(gr)
end(gr)
width(gr)

# the mcols() accessor extracts the metadata as a DataFrame
mcols(gr)

mcols(gr)$score

score(gr)

# add the sequence information (ensure data integrity and prevents accidental mixing of ranges from 
# incompatible context)
seqinfo(gr) <- Seqinfo(genome = 'hg38')
seqinfo(gr)

names(gr)
length(gr)

# subsetting Granges objects
gr[2:3]
gr[2:3,"GC"]
subset(gr,strand == "+" & score > 5, select = score)

# assign elements to the GRanges object
grMod <- gr
grMod[2] <- gr[1]
head(grMod, n = 3)

# methods to repeat, reverse, or select specific portions of GRanges objects
rep(gr[2],times = 3)
rev(gr)
head(gr,n=3)
tail(gr,n=2)
window(gr,start = 2, end = 4)
gr[IRanges(start = c(2,7), end = c(3,9))]


# splitteing and combining GRanges objects
sp <- split(gr,rep(1:2,each = 5))
sp

split(gr,~strand)
c(sp[[1]],sp[[2]])

stack(sp,index.var = "group") #stacks the elements of a GRangeList into a single GrRanges and adds a column
                              #indicating the origin of each element


# aggregating GRanges objects
aggregate(gr,score ~ strand, mean)
aggregate(gr,~strand, n_score = lengths(score),mean_score = mean(score)) # need to call lengths(score) 
                                                                         # instead of length(score)

# basic interval operations for GRanges objects
# flank function can be used to recover regions flanking the set of ranges represented by the GRanges object
g <- gr[1:3]
g <- append(g, gr[10]) 
flank(g, 10) #get a GRange object containing the ranges that include the 10 base 
             #upstream according to the direction of "transcription"
flank(g,10,start = FALSE)

# generate a region starting 2000bp upstream and 200bp downstream of the TSS
promoters(g)

# ignore strand/transcription and assume the orientation of left to right 
flank(unstrand(g), 10)

# move the ranges by a specific number of base pairs
shift(g,5)

# set a specific width, by default fixing the "transcription" start
resize(g,30)

# inter-range functions
# merge overlapping and adjacent ranges to produce a minimal set of ranges representing
# the regions covered by the original set
reduce(gr)
reduce(gr,ignore.strand = TRUE)

gaps(g)

# break up the ranges so that they do not overlap but still cover the same regions
disjoin(g)

# count how many ranges overlap each position in the sequence universe of a GRange object
cov <- coverage(g)
cov_gr <- GRanges(cov)
cov <- coverage(cov_gr, weight = "score")

# compact representation of width 1 ranges
GPos(cov[1:3])


# interval set operations for GRanges objects
g2 <- head(gr,n=2)
union(g,g2)
intersect(g,g2)
setdiff(g,g2)

# element-wise operation
g3 <- g[1:2]
ranges(g3[1]) <- IRanges(start = 105, end = 112)
punion(g2,g3)
pintersect(g2,g3)
psetdiff(g2,g3)

## finding overlaps between GRanges objects
# generate three random data.frame objects
set.seed(66+105+111+99+49+56)

pos <- sample(1:200, size = 30L)
size <- 10L
end <- size + pos - 1L
chrom <- sample(paste0("chr", 1:3), size = 30L, replace = TRUE)
query_df <- data.frame(chrom = chrom, 
                       start = pos,
                       end = end)
query_dfs <- split(query_df, 1:3)
q1 <- rename(query_dfs[[1L]], start = "pos")
q2 <- rename(query_dfs[[2L]], chrom = "ch", start = "st")
q3 <- rename(query_dfs[[3L]], end = "last")

q1 <- makeGRangesFromDataFrame(q1, start.field = "pos")
q2 <- makeGRangesFromDataFrame(q2, seqnames.field = "ch",
                               start.field = "st")
q3 <- makeGRangesFromDataFrame(q3, end.field = "last")
query <- mstack(q1, q2, q3, .index.var="replicate")
sort(query, by = ~ start)

# extracts the elements in the query (the first argument) that overlap at least one element in the subject 
# (the second)
subject <- gr
subsetByOverlaps(query, subject, ignore.strand=TRUE)

hits <- findOverlaps(query, subject, ignore.strand=TRUE)

joined <- query[queryHits(hits)]
joined$score <- subject$score[subjectHits(hits)]

ranges(joined) <- ranges(pintersect(joined, subject[subjectHits(hits)]))

# group the subject hits by query hits
hitsByQuery <- as(hits,"List")

# hitsByQuery is an IntegerList, it has many methods for efficient aggregation
counts <- countQueryHits(hits)
counts <- countOverlaps(query,subject,ignore.strand = TRUE)
unname(counts)

# annotate each query with the maximum score among the subject hits
query$maxScore <- max(extractList(subject$score,hitsByQuery))
subset(query,maxScore > 0)

hits <- findOverlaps(query, subject, select="first", ignore.strand=TRUE)
hits <- findOverlaps(query, subject, select="arbitrary", ignore.strand=TRUE)
hits

# Example: exploring BigWig files from AnnotationHub
library(AnnotationHub)
ah <- AnnotationHub()
roadmap_hub <- query(ah,"EpigenomeRoadMap")
metadata <- query(ah,"Metadata")[[1L]]
head(metadata)

# find out the name of the sample corresponding to primary memory T-cells
primary_tcells <- subset(metadata,
                         ANATOMY == "BLOOD" & TYPE == "PrimaryCell" &
                          EDACC_NAME == "CD8_Memory_Primary_Cells")$EID
# extract the sample ID corresponding to the filter
primary_tcells <- as.character(primary_tcells)

# take roadmap hub and query it based on other conditions
methylation_files <- query(roadmap_hub,
                           c("BigWig",primary_tcells,"H3K4ME[1-3]",
                             "pval.signal"))
methylation_files

bw_files <- lapply(methylation_files[1:2],`[[`,1L)

# extract the genome information from the first BigWig file and filter to get the range for chromosome 10
chr10_ranges <- Seqinfo(genome="hg19")["chr10"]

# read the BigWig file only extracting scores if they overlap chromosome 10
library(rtracklayer)
chr10_scores <- lapply(bw_files, import, which = chr10_ranges,
                       as = "RleList") 
chr10_scores[[1]]$chr10

islands <- lapply(chr10_scores, slice, lower=1L)
summits <- lapply(islands, viewRangeMaxs)

summits <- lapply(lapply(summits, `+`, 50L), reduce)

summits_grs <- lapply(summits, GRanges)
score_grs <- mapply(function(scores, summits) {
  summits$score <- scores[summits]
  seqlengths(summits) <- lengths(scores)
  summits
}, chr10_scores, summits_grs)
score_gr <- stack(GenomicRangesList(score_grs), index.var="signal_type")

score_gr$score_max <- max(score_gr$score)
chr10_max_score_region <- aggregate(score_gr, score_max ~ signal_type, max)

# Example: coverage analysis of BAM files
library(tools)
bams <- list_files_with_exts(system.file("extdata", package = "airway"), "bam")
names(bams) <- sub("_[^_]+$", "", basename(bams))
library(Rsamtools)
bams <- BamFileList(bams)

# compute the coverage of the alignments over all contigs in the BAM
first_bam <- bams[[1L]]
first_bam_cvg <- coverage(first_bam)

head(table(first_bam_cvg)[1L,])

# read the BAM file into a GAlignments object using readGAlignments() and extract the ranges, chopping by introns, using grglist()
library(GenomicAlignments)
reads <- grglist(readGAlignments(first_bam))

reads[lengths(reads) >= 2L]

# count how many reads overlap each gene
# get the transcript structures as a GrangeList from Ensembl
reads[lengths(reads) >= 2L]
library(EnsDb.Hsapiens.v75)
tx <- exonsBy(EnsDb.Hsapiens.v75, "gene")

reads <- keepStandardChromosomes(reads)
counts <- countOverlaps(tx, reads, ignore.strand=TRUE)
head(counts[counts > 0])

airway <- summarizeOverlaps(features=tx, reads=bams,
                            mode="Union", singleEnd=FALSE,
                            ignore.strand=TRUE, fragments=TRUE)