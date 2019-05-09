## TxDb exercises

# How many transcripts does PPARG have, according to UCSC?
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txs <- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene)
txs$`5468`

# Does Ensembl agree?
library(EnsDb.Hsapiens.v79)
txs2 <- transcriptsBy(EnsDb.Hsapiens.v79)
txs2$ENSG00000132170

# How many genes are between 2858473 and 3271812 on chr2 in the hg19 genome?
library(GenomicRanges)
gns <- GRanges(seqnames = "chr2", ranges = IRanges(2858473,3271812))


## OrganismDb exercises

# Get all the GO terms for BRCA1
library(Homo.sapiens)
select(Homo.sapiens,keys = '672',columns = c("GO","GOALL","GOID"),keytype = "ENTREZID")

# What gene does the UCSC transcript ID uc002fai.3 map to?
select(Homo.sapiens,keys = 'uc002fai.3',columns = c('ENTREZID','GENENAME','GENEID'),keytype = 'TXNAME')

# Get all the transcripts from the hg19 genome build, along with their Ensembl gene ID, 
# UCSC transcript ID and gene symbol
transcripts(Homo.sapiens,columns = c("GENEID","TXID","SYMBOL"))


## Organism.dplyr exercises
# How many supported organisms are implemented in Organism.dplyr?
library(Organism.dplyr)
supportedOrganisms()

# Display the ensembl Id and genename description for symbol “NAT2”.
src <- src_ucsc("human")
select_tbl(src,'NAT2',c('ensembl','genename'),'symbol')


## BSgenome 
# Get the sequences for all transcripts of the TP53 gene
library(BSgenome.Hsapiens.UCSC.hg19)
gn <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
getSeq(Hsapiens,gn["7157"])

