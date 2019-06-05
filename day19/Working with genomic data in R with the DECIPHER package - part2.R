## predict homology between genes that are linked by syntenic hits

# if our synteny object was built from 5 genomes, that object is a 5*5 matrix. NucleotideOverlap accesses the upper triangle
# of that object to build a 5*5 matrix where each position is built from data in the analogous position from the Synteny
# object
MatrixObject <- NucleotideOverlap(SyntenyObject = SyntenyObject,
                                  GeneCalls = GeneCalls,
                                  Verbose = TRUE)

# Catalog will return every agglomerated set of pairs, of every size possible for the given set of genomes
Homologs <- Catalog(MatrixObject,
                    Verbose = TRUE)

# visualize this object as a histogram of the size of these agglomerations, by the number of pairs included in each agglomerated
# group
hist(sapply(Homologs,
            function(x) nrow(x)),
     main = "Size of Agglomerations",
     ylab = "Number of Sets",
     xlab = "Number of Gene Pairs",
     breaks = 27L)

# query this list for any other presence absence pattern 
MaxRows <- max(sapply(Homologs, function(x) nrow(x)))

CoreSet <- which(sapply(Homologs, function(x) nrow(x)) == MaxRows)

# CoreAligner function collects the genes in these sets, from their respective genomes, aligns the respective sets, and then
# concatenates the alignments.
CoreGenome <- CoreAligner(Homologs[CoreSet],
                          PATH = DBPath,
                          GeneCalls = GeneCalls,
                          Verbose = TRUE)

CoreDist <- DistanceMatrix(myXStringSet = CoreGenome,
                           verbose = FALSE,
                           correction = "Jukes-Cantor")

CoreDend <- IdClusters(myDistMatrix = CoreDist,
                       myXStringSet = CoreGenome,
                       method = "NJ",
                       verbose = FALSE,
                       showPlot = TRUE,
                       type = "dendrogram")

unlink(DBPath)

# LogicPan function allows us to utilize the Homologs object, and the list of GeneCall dataframes, to create a presence absence matrix
# of the pan genome
PanGenomeMatrix <- LogicalPan(HomologList = Homologs,
                              GeneCalls = GeneCalls,
                              Verbose = TRUE,
                              Plot = FALSE)

# visualize this matrix
image(t(PanGenomeMatrix), col = c("white", "blue"),
      main = "Presence Absence")

# we can create a dendrogram from this matrix
PanGenome <- dist(PanGenomeMatrix,
                  method = "binary")

PanDend <- IdClusters(myDistMatrix = PanGenome,
                      method = "NJ",
                      type = "dendrogram",
                      showPlot = TRUE,
                      verbose = FALSE)

# we can create a simple tangleogram from these two phylogenetic trees
tf1 <- tempfile()
tf2 <- tempfile()

WriteDendrogram(x = PanDend,
                file = tf1)

WriteDendrogram(x = CoreDend,
                file = tf2)

unlink(tf1)
unlink(tf2)
layout(matrix(1:2, nrow = 1L))
p <- par(mar = c(5, 2, 1, 2))
plot(CoreDend, horiz = TRUE, leftlab = "none")
par(mar = c(5, 2, 1, 2))
plot(PanDend, horiz = TRUE, xlim = c(0, attr(PanDend, "height")), leaflab = "none")
C.Ord <- unlist(CoreDend)
P.Ord <- unlist(PanDend)
segments(-0.10, seq_along(C.Ord), -0.01, match(C.Ord, P.Ord), xpd = NA, col = "blue")

par(p)

# collect the annotations for all of the sets of homologs that have been predicted, and see how much these annotations agree, and if
# any of these annotations provide evidence for our proposition
CoreGenes <- matrix(data = NA_integer_,
                    ncol = length(GeneCalls),
                    nrow = length(CoreSet))
CoreAnnotations <- matrix(data = NA_character_,
                          ncol = length(GeneCalls),
                          nrow = length(CoreSet))
for (i in seq_len(ncol(CoreGenes))) {
  for (j in seq_len(nrow(CoreGenes))) {
    CoreGenes[j, i] <- unique(Homologs[CoreSet[j]][[1]][, i])[!is.na(unique(Homologs[CoreSet[j]][[1]][, i]))]
    CoreAnnotations[j, i] <- GeneCalls[[i]][CoreGenes[j, i], "Annotation"]
  }
}

CoreAnnotations <- t(CoreAnnotations)
CoreAnnotations <- data.frame(CoreAnnotations,
                              stringsAsFactors = FALSE)

CoreAnnotations[, c(1L, 3L)]

CoreAnnotations[, 2L]

CoreAnnotations[, c(17L, 31L)]