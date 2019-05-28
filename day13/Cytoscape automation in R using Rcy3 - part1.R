library(RCy3)

installation_responses <- c()
cyto_app_toinstall <- c("clustermaker2", "enrichmentmap", "autoannotate", "wordcloud", "stringapp", "aMatReader")

cytoscape_version <- unlist(strsplit( cytoscapeVersionInfo()['cytoscapeVersion'],split = "\\."))
if(length(cytoscape_version) == 3 && as.numeric(cytoscape_version[1]>=3) 
   && as.numeric(cytoscape_version[2]>=7)){
  for(i in 1:length(cyto_app_toinstall)){
    #check to see if the app is installed.  Only install it if it hasn't been installed
    if(!grep(commandsGET(paste("apps status app=\"", cyto_app_toinstall[i],"\"", sep="")), 
             pattern = "status: Installed")){
      installation_response <-commandsGET(paste("apps install app=\"", 
                                                cyto_app_toinstall[i],"\"", sep=""))
      installation_responses <- c(installation_responses,installation_response)
    } else{
      installation_responses <- c(installation_responses,"already installed")
    }
  }
  installation_summary <- data.frame(name = cyto_app_toinstall, 
                                     status = installation_responses)
  
  knitr::kable(list(installation_summary),
               booktabs = TRUE, caption = 'A Summary of automated app installation'
  )
}

# Make sure that Cytoscape is open 
cytoscapePing()

cytoscapeVersionInfo()

# available functions, commands and arguments
help(package = RCy3)
cyrestAPI()
commandsAPI()
commandsHelp('help')

### Cytoscape Basics

# create a Cytoscape network from some basic R objects
nodes <- data.frame(id=c("node 0","node 1", "node 2", "node 3"),
                    group=c("A","A","B","B"),
                    score=as.integer(c(20,10,15,5)),
                    stringsAsFactors = FALSE)
edges <- data.frame(source=c("node 0","node 0","node 0","node 2"),
                    target=c("node 1","node 2","node 3","node 3"),
                    interaction=c("inhibits","interactis","activates","interacts"),
                    weight=c(5.1,3.0,5.2,9.9),
                    stringsAsFactors = FALSE)
nodes
edges
createNetworkFromDataFrames(nodes,edges, title = "my first network",collection = "DataFrame Example")

# get an image of the resulting network and include it in our current analysis 
initial_network_png_file_name <- file.path(getwd(),"my_first_network.png")
if(file.exists(initial_network_png_file_name)){
  #cytoscape hangs waiting for user response if file already exists.  Remove it first
  file.remove(initial_network_png_file_name)
} 
exportImage(initial_network_png_file_name, type = "png")

# example data set
RNASeq_expression_matrix <- read.table( 
  file.path(getwd(),"TCGA_OV_RNAseq_expression.txt"),  
  header = TRUE, sep = "\t", quote="\"", stringsAsFactors = FALSE)

RNASeq_gene_scores <- read.table( 
  file.path(getwd(),"TCGA_OV_RNAseq_All_edgeR_scores.txt"),  
  header = TRUE, sep = "\t", quote="\"", stringsAsFactors = FALSE)                                      

## Use case 1 - how are my top genes related?
top_mesenchymal_genes <- RNASeq_gene_scores[which(RNASeq_gene_scores$FDR.mesen < 0.05 & RNASeq_gene_scores$logFC.mesen > 2),]
head(top_mensenchymal_genes)

# query the String Database to get all interactions found for our set of top Mesenchymal genes
commandsHelp("help string")
commandsHelp("help string protein query")

mesen_string_interaction_cmd <- paste('string protein query taxonID=9606 limit=150 cutoff=0.9 query="',paste(top_mesenchymal_genes$Name, collapse=","),'"',sep="")
commandsGET(mesen_string_interaction_cmd)

initial_string_network_png_file_name <- file.path(getwd(), "initial_string_network.png")
if(file.exists(initial_string_network_png_file_name)){
  #cytoscape hangs waiting for user response if file already exists.  Remove it first
  response <- file.remove(initial_string_network_png_file_name)
} 

response <- exportImage(initial_string_network_png_file_name, type = "png")

layoutNetwork('force-directed')

getLayoutNames()
getLayoutPropertyNames(layout.name='force-directed')
layoutNetwork('force-directed defaultSpringCoefficient=0.0000008 defaultSpringLength=70')


