## Workflow for multi-omics analysis with MultiAssayExperiment

library(MultiAssayExperiment)
library(GenomicRanges)
library(RaggedExperiment)
library(curatedTCGAData)
library(GenomicDataCommons)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(TCGAutils)
library(UpSetR)
library(mirbase.db)
library(AnnotationFilter)
library(EnsDb.Hsapiens.v86)
library(survival)
library(survminer)
library(pheatmap)

# Overview of key data classes

# 1. (Ranged) SummarizedExperiment

## SummarizedExperiment can store multiple experimental data matrices of identical dimensions, with associated metadata on the rows/
## genes/transcripts/other measurements, column/sample phenotype or clinical data, and the overall experiment.
## RanngedSummarizedExperiment is a derivative class, and SingleCellExperiment extends RangedSummarizedExperiment, which in turn extends
## SummarizedExperiment.

library(SingleCellExperiment)
extends("SingleCellExperiment")

# 2. RaggedExperiment

## RaggedExperiment is a flexible data representation for segmented copy number, somatic mutations such as represented 
## in .vcf files, and other ragged array schema for genomic location data. Like the GRangesList class from GenomicRanges, 
## RaggedExperiment can be used to represent differing genomic ranges on each of a set of samples.

showClass("RaggedExperiment")

# 3. MultiAssayExperiment

## Integrative container for coordinating multi-omics experiment data on a set of biological specimens.

################################################

# Working with RaggedExperiment

# 1. constructing a RaggedExperiment object

# we starte with a toy example of two GRanges objects, providing ranges on 2 chromosomes in two samples
sample1 <- GRanges(c(A = "chr1:1-10:-", B = "chr1:8-14:+", C = "chr1:15-18:+"),
                   score = 3:5, type = c("germline", "somatic", "germline"))
sample2 <- GRanges(c(D = "chr1:1-10:-", E = "chr1:11-18:+"),
                   score = 11:12, type = c("germline", "somatic"))
colDat <- DataFrame(id = 1:2, status = factor(c("control", "case")))

# RaggedExperiment can be constructed from individual Granges:
(ragexp <- RaggedExperiment(
  sample1 = sample1,
  sample2 = sample2,
  colData = colDat))

# or from a GRangeList
grl <- GRangesList(sample1 = sample1, sample2 = sample2)
ragexp2 <- RaggedExperiment(grl, colData = colDat)

# the original ranges are represented as the rowRanges of the RaggedExperiment
rowRanges(ragexp)
