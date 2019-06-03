## Data import and creating pipelines

# Worked example: exploring BigWig files from AnnotationHub

# extracting data from AnnotationHub

library(AnnotationHub)
library(magrittr)
ah <- AnnotationHub() # first we construct a hub that contains all references to the EpigenomeRoadMap data and extract the
                      # metadata as a data.frame

roadmap_hub <- ah %>% 
  query("EpigenomeRoadMap")

metadata <- ah %>%
  query("Metadata") %>%
  extract2(names(.))

head(metadata)

primary_tcells <-  metadata %>% 
  filter(ANATOMY == "BLOOD") %>% 
  filter(TYPE == "PrimaryCell") %>% 
  filter(EDACC_NAME == "CD8_Memory_Primary_Cells") %>% 
  extract2("EID") %>% 
  as.character()
primary_tcells

methylation_files <- roadmap_hub %>% 
  query("BigWig") %>% 
  query(primary_tcells) %>% 
  query("H3K4ME[1-3]") %>% 
  query("pval.signal")

methylation_files

bw_files <- lapply(c("AH33454", "AH33455"), function(id) ah[[id]]) 

names(bw_files) <- c("HK34ME1", "HK34ME3")

### reading BigWig files
chr10_ranges <- bw_files %>% 
  extract2(1L) %>% 
  get_genome_info() %>%
  filter(seqnames == "chr10")

# a function to read our bigwig files
read_chr10_scores <- function(file) {
  read_bigwig(file, overlap_ranges = chr10_ranges) %>% 
    set_genome_info(genome = "hg19")
}
# apply the function to each file
chr10_scores <- lapply(bw_files, read_chr10_scores) 
# bind the ranges to a single GRAnges object
# and add a column to identify the ranges by signal type
chr10_scores <- bind_ranges(chr10_scores, .id = "signal_type")
chr10_scores

# since we want to call peaks over each signal type, we will create a grouped GRanges object
chr10_scores_by_signal <- chr10_scores %>% 
  group_by(signal_type)

# filter to find the coordinates of the peak containing the maximum score for each signal
chr10_max_score_region <- chr10_scores_by_signal %>%
  filter(score == max(score)) %>% 
  ungroup() %>% 
  anchor_center() %>%
  mutate(width = 5000)

peak_region <- chr10_scores %>%
  join_overlap_inner(chr10_max_score_region) %>% 
  filter(signal_type.x == signal_type.y) 


### coverage analysis of BAM files

# first gather all the BAM files available to use in airway
bfs <- system.file("extdata", package = "airway") %>% 
  dir(pattern = ".bam", 
      full.names = TRUE)

names(bfs) <- bfs %>% 
  basename() %>% 
  sub("_[^_]+$", "", .)

# compute the coverage of the alignments over all contigs in the BAM
first_bam_cvg <- bfs %>% 
  extract2(1) %>% 
  compute_coverage()
first_bam_cvg

first_bam_cvg %>% 
  group_by(seqnames, score) %>% 
  summarise(n_bases_covered = sum(width))

split_bam <- bfs %>% 
  extract2(1) %>% 
  read_bam(index = NULL) 

split_bam

split_bam %>% 
  select(flag)

split_bam <- split_bam %>% 
  chop_by_introns()
split_bam

split_bam %>% 
  filter(n() >= 2)
