suppressMessages(library(DECIPHER))
library(FindHomology)
library(stringr)

# Genomes, and DECIPHER databases
data("GeneCallAdds") # a character vector of ftp addresses for genomic GFF files
data("GenomeAdds")  # a character vector of ftp addresses for genomic fasta files
data("GenomeIDs") # a vector of abbreviated names for the genomes 

GeneCalls <- FindHomology::GFFParser(GFFAddress = GeneCallAdds,
                                     Verbose = TRUE)
# DECIPHER and Databases
# we can access our databases through either a database connection, or simply a filepath to the database
# in the code below, we will use a database connection
DBPath <- tempfile()

DBConn <- dbConnect(SQLite(),
                    DBPath)

for (i in seq_along(GenomeAdds)){
  Seqs2DB(seqs = GenomeAdds[i],
          type = "FASTA",
          dbFile = DBConn,
          identifier = as.character(i),
          tblName = "Seqs",
          verbose = FALSE)
}

# identifying details of the database can be viewed in the default browser
BrowseDB(DBPath)

# Comparison of genomes
# determine Synteny
SyntenyObject <-  FindSynteny(dbFile = DBPath,
                              verbose = TRUE)
plot(SyntenyObject[2:3, 2:3])

plot(SyntenyObject[2:3, 2:3],
     "frequency")

plot(SyntenyObject[2:3, 2:3],
     "neighbor")

plot(SyntenyObject,
     "neighbor")

pairs(SyntenyObject[2:3, 2:3], labels = GenomeIDs[2:3])

pairs(SyntenyObject, labels = GenomeIDs)
