## Working with MultiAssayExperiment
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
miniACC[c("MAPK14", "IGFBP2"), , ]
miniACC[, miniACC$pathologic_stage == "stage iv"]
miniACC[, , "RNASeq2GeneNorm"]
