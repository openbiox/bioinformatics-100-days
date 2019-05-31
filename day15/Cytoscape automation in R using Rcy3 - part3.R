## Use case 3 - Functional Enrichment of Omics set

suppressPackageStartupMessages(library(EnrichmentBrowser))

# download the latest pathway definition file from the Baderlab download site
tryCatch(expr = { suppressPackageStartupMessages(library("RCurl"))}, 
         error = function(e) {  install.packages("RCurl")}, 
         finally = library("RCurl"))

gmt_url = "http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/"


#list all the files on the server
filenames = getURL(gmt_url)
tc = textConnection(filenames)
contents = suppressWarnings(readLines(tc))
close(tc)

#get the gmt that has all the pathways and does not include terms inferred from electronic annotations(IEA)
#start with gmt file that has pathways only
rx = gregexpr("(?<=<a href=\")(.*.GOBP_AllPathways_no_GO_iea.*.)(.gmt)(?=\">)",
              contents, perl = TRUE)
gmt_file = unlist(regmatches(contents, rx))

dest_gmt_file <- file.path(getwd(), gmt_file)

download.file(
  paste(gmt_url,gmt_file,sep=""),
  destfile=dest_gmt_file
)

# load in the gmt file
baderlab.gs <- getGenesets(dest_gmt_file)
baderlab.gs

# create the dataset required by EnrichmentBrowser tools

# create the expression file - A tab separated text file containing expression values. Columns = sample/subjects; row =
# features/probes/genes; No headers, row or column names
expr <- RNASeq_expression

sumexpr_filename <- file.path(getwd(), "SummarizeExperiment_expression.txt")
write.table(expr, file = sumexpr_filename, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

rowData <- RNASeq_gene_scores[,grep(colnames(RNASeq_gene_scores), pattern = "mesen")]
rowData <- cbind(RNASeq_gene_scores$Name, rowData)
colnames(rowData)[2] <- "FC"
colnames(rowData)[6] <- "ADJ.PVAL"

sumexpr_rdat_filename <- file.path(getwd(), "SummarizeExperiment_rdat.txt")
write.table(rowData[,1], file = sumexpr_rdat_filename, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

# load in the data classification data
# A tab separated text file containing annotation information for the samples in either *two or three* columns. No headers,
# row or column names. The number of rows/samples in this file should match the number of columns/samples of the expression
# matrix. The first column is reserved for the sample IDs; The second columns is reserved for a *BINARY* group assignment.
# Use '0' and '1' for unaffected (controls) and affected (cases) sample class
classDefinitions_RNASeq <- read.table(
  file.path(getwd(), "RNASeq_classdefinitions.txt"), header = TRUE, sep = "\t", quote = "\"", stringsAsFactors = FALSE)

colData <- data.frame(Sample = colnames(RNASeq_expression),
                      GROUP = classDefinitions_RNASeq$SUBTYPE,
                      stringsAsFactors = FALSE)
rownames(colData) <- colnames(RNASeq_expression)
colData$GROUP[which(colData$GROUP != "Mesenchymal")] <- 0
colData$GROUP[which(colData$GROUP == "Mesenchymal")] <- 1

sumexpr_cdat_filename <- file.path(getwd(), "SummarizeExperiment_cdat.txt")
write.table(colData, file = sumexpr_cdat_filename, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

# create the Summarize Experiment object
se_OV <- EnrichmentBrowser::readSE(assay.file = sumexpr_filename, cdat.file = sumexpr_cdat_filename, rdat.file = sumexpr_rdat_filename)

# put our precomputed p-values and fold change values into the Summarized Experiment object so we can use our rankings for the analysis
# set the Summarized Experiment to our computed p-values and FC
rowData(se_OV) <- rowData

# run basic Over representation analysis using our ranked genes and our gene set file downloaded from the Baderlab genesets
ora.all <- sbea(method = "ora", se = se_OV, gs = baderlab.gs, perm = 0, alpha = 0.05)
gsRanking(ora.all)

# take the enrichment results and create a generic enrichment map input file so we can create an Enrichment map.

# manually adjust p-values
ora.all$res.tbl <- cbind(ora.all$res.tbl, p.adjust(ora.all$res.tbl$PVAL, "BH"))
colnames(ora.all$res.tbl)[ncol(ora.all$res.tbl)] <- "Q.VALUE"

# create a generic enrichment map file
em_results_mesen <- data.frame(name = ora.all$res.tbl$GENE.SET, descr = ora.all$res.tbl$GENE.SET, 
                               pvalue = ora.all$res.tbl$PVAL, qvalue = ora.all$res.tbl$Q.VALUE, stringsAsFactors = FALSE)

# write out the ora results
em_results_mesen_filename <- file.path(getwd(), "mesen_ora_enr_results.txt")

write.table(em_results_mesen,em_results_mesen_filename,col.name = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)

# create an enrichment map with the returned ORA results
em_command = paste('enrichmentmap build analysisType="generic" ', "gmtFile=", dest_gmt_file,
                   'pvalue=',"0.05", 'qvalue=',"0.05",
                   'similaritycutoff=',"0.25",
                   'coeffecients=',"JACCARD",
                   'enrichmentsDataset1=',em_results_mesen_filename ,
                   sep=" ")

em_mesen_network_suid <- commandsRun(em_command)

renameNetwork("Mesenchymal_ORA_enrichmentmap", network = as.numeric(em_mesen_network_suid))

mesenem_png_file_name <- file.path(getwd(), "mesenem.png")
exportImage(mesenem_png_file_name, type = "png")

# annotate the enrichment map to get the general themes that are found in the enrichment results

#get the column from the nodetable
nodetable_colnames <- getTableColumnNames(table="node",  network =  as.numeric(em_mesen_network_suid))

descr_attrib <- nodetable_colnames[grep(nodetable_colnames, pattern = "_GS_DESCR")]

#Autoannotate the network
autoannotate_url <- paste("autoannotate annotate-clusterBoosted labelColumn=", descr_attrib," maxWords=3 ", sep="")
current_name <-commandsGET(autoannotate_url)

#create a summary network 
commandsGET("autoannotate summary network='current'")

#change the network name
summary_network_suid <- getNetworkSuid()

renameNetwork(title = "Mesen_ORA_summary_netowrk",
              network = as.numeric(summary_network_suid))

#get the summary node names
summary_nodes <- getTableColumns(table="node",columns=c("name"))

#clear bypass style the summary network has
clearNodePropertyBypass(node.names = summary_nodes$name,visual.property = "NODE_SIZE")

#relayout network
layoutNetwork('cose',
              network = as.numeric(summary_network_suid))

mesenem_summary_png_file_name <- file.path(getwd(),"mesenem_summary_network.png")

exportImage(mesenem_summary_png_file_name, type = "png")
