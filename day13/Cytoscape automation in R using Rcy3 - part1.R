library(RCy3)

# if you use Cytoscape 3.7 or higher then apps can be installed directly from R
# if you use Cytoscape older than 3.7, you need to install these apps manually from app store
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
commandsHelp("help")

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

# idea:
# query the String Database to get all interactions found for our set of top Mesenchymal genes
# StrinApp in cytoscape is a protein-protein and protein-chemical database that imports data from String(Szklarczyk et al. 2016)
# into s unified, queriable database

# get a subset of genes of interest from our scored data
top_mesenchymal_genes <- RNASeq_gene_scores[which(RNASeq_gene_scores$FDR.mesen < 0.05 & RNASeq_gene_scores$logFC.mesen > 2),]
head(top_mesenchymal_genes)

# we are going to query the String Database to get all interactions found for our set of top Mesenchymal genes
# to see the parameters required by the string function or to find the right string function we can use
commandsHelp("help string")
commandsHelp("help string protein query")

mesen_string_interaction_cmd <- paste('string protein query taxonID=9606 limit=150 cutoff=0.9 query="',paste(top_mesenchymal_genes$Name, collapse=","),'"',sep="")
mesen_string_interaction_cmd
commandsGET(mesen_string_interaction_cmd)

# get a screenshot of the initial network
initial_string_network_png_file_name <- file.path(getwd(), "initial_string_network.png")
if(file.exists(initial_string_network_png_file_name)){
  #cytoscape hangs waiting for user response if file already exists.  Remove it first
  response <- file.remove(initial_string_network_png_file_name)
} 

response <- exportImage(initial_string_network_png_file_name, type = "png")

# layout the network
layoutNetwork('force-directed')

# check what other layout algorithms are available to try out
getLayoutNames()

# get the parameters for a specific layout
getLayoutPropertyNames(layout.name='force-directed')

# re-layout the network using the force directed layout but specify some of the parameters
layoutNetwork('force-directed defaultSpringCoefficient=0.0000008 defaultSpringLength=70')

# get a screenshot of the re-laid out netword
relayout_string_network_png_file_name <- file.path(getwd(),"relayout_string_network.png")
if(file.exists(relayout_string_network_png_file_name)){
  #cytoscape hangs waiting for user response if file already exists.  Remove it first
  response<- file.remove(relayout_string_network_png_file_name)
} 
response <- exportImage(relayout_string_network_png_file_name, type = "png")
getTableColumnNames('node')

node_attribute_table_topmesen <- getTableColumns(table="node")
head(node_attribute_table_topmesen[,3:7])

## overlay our expression data on the String network
# check what identifiers type is used by String by pulling in the column names of the node attribute table
getTableColumnNames('node')

# pull in the entire node attribute table
node_attribute_table_topmesen <- getTableColumns(table = "node")
head(node_attribute_table_topmesen[,3:7])

# the column "display name" contains HGNC names which are also found in our Ovarian Cancer dataset
# to import our expression data we will match our dataset to the "display name" node attribute
?loadTableData
loadTableData(RNASeq_gene_scores, table.key.column = "display name", data.key.column = "Name")

# modify the visual style 
# start with a default style
style.name = "MesenchymalStyle"
default.list <- list(NODE_SHAPE = "eclipse",
                     NODE_SIZE = 60,
                     NODE_FILL_COLOR = "#AAAAAA",
                     EDGE_TRANSPARENCY = 120)
node.label.map <- mapVisualProperty('node label','display name','p')
createVisualStyle(style.name,default.list,list(node.label.map))
setVisualStyle(style.name = style.name)

# grab the column data from Cytoscape and pull out the min and max to define our data mapping range of values
min.mesen.logfc = min(RNASeq_gene_scores$logFC.mesen,na.rm = TRUE)
max.mesen.logfc = max(RNASeq_gene_scores$logFC.mesen,na.rm = TRUE)
data.values = c(min.mesen.logfc,0,max.mesen.logfc)

# use the RColorBrewer package to help us pick good colors to pair with our data values
library(RColorBrewer)
display.brewer.all(length(data.values),colorblindFriendly = TRUE,type = "div")
node.colors <- c(rev(brewer.pal(length(data.values),"RdBu")))

# map the colos to our data value and update our visual style
setNodeColorMapping("logFC.mesen",data.values, node.colors,style.name=style.name)

# String includes our query proteins as well as other proteins that associate with our query proteins
# not all of the proteins in this network are our tp hits
# we want to visualize the proteins that  are our top Mesenchymal hits

# add a different border color or change the node shape for our top hits
getNodeShapes()
setNodeShapeBypass(node.names = top_mesenchymal_genes$Name, new.shapes = "TRIANGLE")

# change the size of the node to be correlated with the Mesenchymal p-value
setNodeSizeMapping(table.column = 'LR.mesen',
                   table.column.values = c(min(RNASeq_gene_scores$LR.mesen),
                                           mean(RNASeq_gene_scores$LR.mesen),
                                           max(RNASeq_gene_scores$LR.mesen)),
                   sizes = c(30,60,150),mapping.type = "c",style.name = style.name)

mesen_string_network_png_file_name <- file.path(getwd(),"mesen_string_network.png")
if(file.exists(mesen_string_network_png_file_name)){
  #cytoscape hangs waiting for user response if file already exists.  Remove it first
  response<- file.remove(mesen_string_network_png_file_name)
} 
response <- exportImage(mesen_string_network_png_file_name, type = "png")
