# Use Case 2 - which genes have similar expression

# we can transform our expression dataset into a correlation matrix
# Using the Cytoscape App, aMatReader, we can transform our adjacency matrix into an interaction network.

# first we filter the correaltion matrix to contain only the strongest connections
RNASeq_expression <- RNASeq_expression_matrix[,3:ncol(RNASeq_expression_matrix)]
RNASeq_expression
rownames(RNASeq_expression) <- RNASeq_expression_matrix$Name 
RNASeq_correlation_matrix <- cor(t(RNASeq_expression), method = "pearson")

# set the diagonal of matrix to zero - eliminate self-correlation
RNASeq_correlation_matrix[row(RNASeq_correlation_matrix) == col(RNASeq_correlation_matrix)] <- 0

# set all correlations that are less than 0.9 to zero
RNASeq_correlation_matrix[which(RNASeq_correlation_matrix<0.90)] <- 0

# get rid of row and columns that have no correlations with the above thresholds
RNASeq_correlation_matrix <- RNASeq_correlation_matrix[which(rowSums(RNASeq_correlation_matrix) != 0),
                                                       which(colSums(RNASeq_correlation_matrix) != 0)]

# write out the correlation file
correlation_filename <- file.path(getwd(),"TCGA_OV_RNAseq_expression_correlation_matrix.txt")
write.table(RNASeq_correlation_matrix,file = correlation_filename, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# use the CyRest call to access the aMatReader functionality
amat_url <- "aMatReader/v1/import"
amat_params = list(files = list(correlation_filename),
                   delimiter = "TAB",
                   undirected = FALSE,
                   ignoreZeros = TRUE, 
                   interactionName = "correlated with",
                   rowNames = FALSE)

response <- cyrestPOST(operation = amat_url, body = amat_params, base.url = "http://localhost:1234")
current_network_id <- response$data["suid"]

#relayout network
layoutNetwork('cose',
              network = as.numeric(current_network_id))

renameNetwork(title = "Coexpression_network_pear0_95",
              network = as.numeric(current_network_id))

# modify the visualization to see where each genes is predominantly expressed
# look at the 4 different p-values associated with each gene and color the nodes with the type associated with the lowest FDR
# load in the scoring data, specify the cancer type where the genes has the lowest FDR values
nodes_in_network <- rownames(RNASeq_correlation_matrix)

# add an additional column to the gene scores table to indicate in which samples the gene is significant
node_class <- vector(length = length(nodes_in_network),mode = "character")
for(i in 1:length(nodes_in_network)){
  current_row <- which(RNASeq_gene_scores$Name == nodes_in_network[i])
  min_pvalue <- min(RNASeq_gene_scores[current_row,
                                       grep(colnames(RNASeq_gene_scores), pattern = "FDR")])
  if(RNASeq_gene_scores$FDR.mesen[current_row] <=min_pvalue){
    node_class[i] <- paste(node_class[i],"mesen",sep = " ")
  }
  if(RNASeq_gene_scores$FDR.diff[current_row] <=min_pvalue){
    node_class[i] <- paste(node_class[i],"diff",sep = " ")
  }
  if(RNASeq_gene_scores$FDR.prolif[current_row] <=min_pvalue){
    node_class[i] <- paste(node_class[i],"prolif",sep = " ")
  }
  if(RNASeq_gene_scores$FDR.immuno[current_row] <=min_pvalue){
    node_class[i] <- paste(node_class[i],"immuno",sep = " ")
  }
}
node_class <- trimws(node_class)
node_class_df <-data.frame(name=nodes_in_network, node_class,stringsAsFactors = FALSE)

head(node_class_df)

# map the new node attribute and the gene scores to the network
loadTableData(RNASeq_gene_scores, table.key.column = "name", data.key.column = "Name")
loadTableData(node_class_df,table.key.column = "name",data.key.column = "name")

# Create a color mapping for the different cancer types
unique_types <- sort(unique(node_class))
coul = brewer.pal(4, "Set1")

coul = colorRampPalette(coul)(length(unique_types))
setNodeColorMapping(table.column = "node_class", table.column.values = unique_types,
                    colors = coul,mapping.type = "d")

correlation_network_png_file_name <- file.path(getwd(),"correaltion_network.png")

if(file.exists(correlation_network_png_file_name)){
  #cytoscape hangs waiting for user response if file already exists.  Remove it first
  file.remove(correlation_network_png_file_name)
} 

exportImage(correlation_network_png_file_name,type = "png")

# cluster the network
setCurrentNetwork(network = getNetworkName(suid = as.numeric(current_network_id)))

clustermaker_url <- paste("cluster mcl network=SUID:", current_network_id, sep = "")
commandsGET(clustermaker_url)

default_node_table <- getTableColumns(table = "node", network = as.numeric(current_network_id))
head(default_node_table)

# perform pathway enrichment on one of the clusters using g:Profiler
tryCatch(expr = {library("gProfileR")},
         error = function(e) {install.packages("gProfileR")}, finally = library("gProfileR"))

runGprofiler <- function(genes,current_organism = "hsapiens", 
                         significant_only = F, set_size_max = 200, 
                         set_size_min = 3, filter_gs_size_min = 5 , exclude_iea = F){
  
  gprofiler_results <- gprofiler(genes ,
                                 significant=significant_only,ordered_query = F,
                                 exclude_iea=exclude_iea,max_set_size = set_size_max,
                                 min_set_size = set_size_min,
                                 correction_method = "fdr",
                                 organism = current_organism,
                                 src_filter = c("GO:BP","REAC"))
  
  #filter results
  gprofiler_results <- gprofiler_results[which(gprofiler_results[,'term.size'] >= 3
                                               & gprofiler_results[,'overlap.size'] >= filter_gs_size_min ),]
  
  # gProfileR returns corrected p-values only.  Set p-value to corrected p-value
  if(dim(gprofiler_results)[1] > 0){
    em_results <- cbind(gprofiler_results[,
                                          c("term.id","term.name","p.value","p.value")], 1,
                        gprofiler_results[,"intersection"])
    colnames(em_results) <- c("Name","Description", "pvalue","qvalue","phenotype","genes")
    
    return(em_results)
  } else {
    return("no gprofiler results for supplied query")
  }
}

# run g:Profiler, it will return a set of pathways and functions that are found to be enriched in our query set of genes
current_cluster <- "1"
selectednodes <- selectNodes(current_cluster,by.col = "__mclCluster")

# create a subnetwork with cluster 1
subnetwork_suid <- createSubnetwork(nodes = "selected")

renameNetwork("Cluster1_Subnetwork", network = as.numeric(subnetwork_suid))

subnetwork_node_table <- getTableColumns(table = "node", network = as.numeric(subnetwork_suid))

em_results <- runGprofiler(subnetwork_node_table$name)

em_results_filename <- file.path(getwd(),paste("gprofile_cluster",current_cluster,"enr_results.txt",sep = "_"))

write.table(em_results,em_results_filename,col.name = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)

head(em_results)

# create an enrichment map with the returned results
# an enrichment map is a different sort of network, nodes represent pathways or functions, edges between these pathways
# or functions represent shared genes or pathway crosstalk

em_command = paste('enrichmentmap build analysisType="generic" ', 
                   'pvalue=',"0.05", 'qvalue=',"0.05",
                   'similaritycutoff=',"0.25",
                   'coeffecients=',"JACCARD",
                   'enrichmentsDataset1=',em_results_filename ,
                   sep=" ")

em_network_suid <- commandsRun(em_command)

renameNetwork("Cluster1_enrichmentmap", network = as.numeric(em_network_suid))

cluster1em_png_file_name <- file.path(getwd(),"cluster1em.png")

if(file.exists(cluster1em_png_file_name)){
  #cytoscape hangs waiting for user response if file already exists.  Remove it first
  file.remove(cluster1em_png_file_name)
} 

exportImage(cluster1em_png_file_name,type = "png")

# annotate the enrichment map to get the general themes that are found in the enrichment results of cluster 1
nodetable_colnames <- getTableColumnNames(table = "node", network = as.numeric(em_network_suid))

descr_attrib <- nodetable_colnames[grep(nodetable_colnames, pattern = "_GS_DESCR")]

autoannotate_url <- paste("autoannotate annotate-clusterBoosted labelColumn=", descr_attrib," maxWords=3 ", sep="")
current_name <-commandsGET(autoannotate_url)

cluster1em_annot_png_file_name <- file.path(getwd(),"cluster1em_annot.png")

if(file.exists(cluster1em_annot_png_file_name)){
  #cytoscape hangs waiting for user response if file already exists.  Remove it first
  file.remove(cluster1em_annot_png_file_name)
}

exportImage(cluster1em_annot_png_file_name, type = "png")
