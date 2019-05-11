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
