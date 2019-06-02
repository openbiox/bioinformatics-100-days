# Fluent genomic data analysis with plyranges

library(BiocManager)
BiocManager::install(c("plyranges","AnnotationHub","airway"))

# since GRanges object is similar to a data.frame, we can use plyranges to construct a GRanges object from a data.frame
library(plyranges, quietly = TRUE)
genes <- data.frame(seqnames = 'VI',
                    start = c(3322, 3030, 1437, 5066, 6426, 836),
                    end = c(3846, 3338, 2615, 5521, 7565, 1363),
                    strand = c("-","-","-","+","+","+"),
                    gene_id = c("YFL064C","YFL065C","YFL066C",
                                "YFL063W","YFL062W","YFL067W"),
                    stringsAsFactors = FALSE)
gr <- as_granges(genes)
gr

# the grammar
# core verbs
set.seed(2019-06-02)
gr2 <- gr %>%
  mutate(gene_type = "ORF",  # mutate() function is used to add columns
         gc_content = runif(n())) %>% # n() operator returns the number of ranges in GRanges object
  filter(width > 400)  # filter() operation returns ranges if the expression evaluates to TRUE

gr2

# multile expressions can be composed together and will be evaluated as &
gr2 %>%
  filter(strand == "+", gc_content > 0.5)

gr2 %>%
  filter(strand == "+" & gc_content > 0.5)

gr2 %>%
  filter(strand == "+" | gc_content > 0.5)

gr2 %>%
  summarise(avg_gc = mean(gc_content),
            n = n())

gr2 %>%
  group_by(strand) %>%
  summarise(avg_gc = mean(gc_content),
            n = n())

by_strand <- gr2 %>% group_by(strand)
by_strand

by_strand %>%
  filter(n() > 2)

by_strand %>%
  mutate(avg_gc_strand = mean(gc_content))

by_strand %>%
  ungroup()

gr2 %>%
  select(gene_id, gene_type)

gr2 %>%
  select(-gc_content)

gr2 %>%
  select(1:2)

# verbs specific to GRanges 
# arithmetic
gr %>%
  mutate(width = width + 1)

gr %>%
  anchor_end() %>%
  mutate(width = width * 2)

gr %>%
  anchor_center() %>%
  mutate(width = width * 2)

# Genomic aggregation
# reduce_ranges() and disjoin_ranges()
# the reduce verb merges overlapping and neighbouring ranges
gr %>% reduce_ranges()

# we could find out whcih genes are overlapping each other by aggregating over the gene_id column and storing the result
# in a List column
gr %>%
  reduce_ranges(gene_id = List(gene_id))

# The disjoin verb takes the union of end points over all ranges, and results in an expanded range
gr %>% 
  disjoin_ranges(gene_id = List(gene_id))

# Overlaps
# suppose we have some additional measurements that obtained from a new experiment on yeast
# these measurements are for three different replicated and represent single nucleotide or insertion deletion intensities
# from an array
set.seed(66+105+111+99+49+56)

pos <- sample(1:10000, size = 100)
size <- sample(1:3, size = 100, replace = TRUE)
rep1 <- data.frame(chr = "VI", 
                   pos = pos,
                   size = size,
                   X = rnorm(100, mean = 2),
                   Y = rnorm(100, mean = 1))

rep2 <- data.frame(chrom = "VI", 
                   st = pos,
                   width = size,
                   X = rnorm(100, mean = 0.5, sd = 3),
                   Y = rnorm(100, sd = 2))

rep3 <- data.frame(chromosome = "VI", 
                   start = pos,
                   width = size,
                   X = rnorm(100, mean = 2, sd = 3),
                   Y = rnorm(100, mean = 4, sd = 0.5))

rep1 <- as_granges(rep1, seqnames = chr, start = pos, width = size)
rep2 <- as_granges(rep2, seqnames = chrom, start = st)
rep3 <- as_granges(rep3, seqnames = chromosome)

intensities <- bind_ranges(rep1, rep2, rep3, .id = "replicate")
arrange(intensities, start)

olap <- filter_by_overlaps(intensities, gr)
olap

olap <- join_overlap_inner(intensities, gr)

join_overlap_left(intensities, gr)

gr %>%
  mutate(gene_length = width) %>%
  join_overlap_intersect(gr, suffix = c(".query",".subject")) %>%
  filter(gene_id.query != gene_id.subject) %>%
  mutate(folap = width / gene_length)
