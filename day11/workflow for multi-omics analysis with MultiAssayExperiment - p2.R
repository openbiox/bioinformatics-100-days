## Working with MultiAssayExperiment - part2
# The MultiAssayExperiment miniACC demo object
data("miniACC")
miniACC

# colData - information biological units
# This slot is a DataFrame describing the characteristics of biological units.
colData(miniACC)[1:4,1:4]
table(miniACC$race)

# ExperimentList - experiment data
experiments(miniACC)

# sampleMap - relationship graph
# a graph representation of the relationship between biological units and experimental results
sampleMap(miniACC)

# metadata
# metadata can be used to keep additional information about patients, assays performed on individuals or on the entire
# cohort, or features such as genes, proteins, and genomic ranges
metadata(miniACC)

##################

# MultiAssayExperiment Subsetting
# single bracket
# multiassayexperiment[i = rownames, j = primary or colnames, k = assay]
miniACC[c("MAPK14", "IGFBP2"), , ]  # return any row named "MAPK14" or "IGFBP2"
miniACC[, miniACC$pathologic_stage == "stage iv"] # keep only patients of pathological stage iv, and all their associated assays
miniACC[, , "RNASeq2GeneNorm"] # keep only the RNA-seq dataset, and only patients for which this assay is available

################

# Subsetting by genomic ranges
# Double bracket [[
# "[[" is a convenience function for extracting a single element of the MultiAssayExperiment ExperimentList
miniACC[[1L]]

################

# Compleye cases
# complete.cases() shows which patients have complete data for all assays
summary(complete.cases(miniACC))

# intersectColumns() will select complete cases and rearrange each ExperimentList element so its columns correspond exactly to rows
# of colData in the same order
accmatched = intersectColumns(miniACC)
colnames(accmatched)

# row names that are common across assays
# intersectRows() keeps only rows that are common to each assay, and aligns them in identical order
accmatched2 <- intersectRows(miniACC[, , c("RNASeq2GeneNorm", "gistict", "Mutations")])
rownames(accmatched2)


#################

# Extraction
# assay and assays
# assay method will extract the first element of the ExperimentList and will return a matrix 
assay(miniACC)
class(assay(miniACC))

# assays method will return a SimpleList of the data with each element being a matrix
assays(miniACC)

# whereas the [[ returned an assay as its original class, assay() and assays() convert the assay data to matrix form

################

# Transformation/ reshaping
# longFormat 
# in long format a single column provides all assay results, with additional optional colData columns whose values are repeated as 
# necessary
longFormat(miniACC[c("TP53","CTNNB1"), , ],
           colDataCols = c("vital_status", "days_to_death"))

# wideFormat
# in wide format, each feature from each assay goes in a separate column, with one row per primary identifier (patient)
wideFormat(miniACC[c("TP35", "CTNNB1"), , ],
           colDataCols = c("vital_status", "days_to_death"))


################
# MultiAssayExperiment class construction and concatenation
# MultiAssayExperiment constructor function
# The MultiAssayExperiment constructor function can take three arguments:
#   1. experiments - An ExperimentList or list of data
#   2. colData - A DataFrame describing the patients (for cell lines, or other biological units)
#   3. sampleMap - A DataFrame of assay, primary, and colname identifiers

MultiAssayExperiment(experiments = experiments(miniACC),
                     colData = colData(miniACC),
                     sampleMap = sampleMap(miniACC),
                     metadata = metadata(miniACC))

# preMultiAssay - Constructor function helper
# perMultiAssay function allows the user to diagnose typical problems when creating a MultiAssayExperiment object

# c- concatenate to MultiAssayExperiment
# c function allows the user to concatenate an additional experiment to an existing MultiAssayExperiment
miniACC2 <- c(miniACC, log2rnaseq = log2(assays(miniACC) $ RNASeq2GeneNorm), mapFrom=1L) # the mapFrom argument allows the user to
                                                                                         # map from a particular experiment provided
                                                                                         # that the order of the colnames is in the same
experiments(miniACC2)



