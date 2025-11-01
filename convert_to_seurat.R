library(Seurat)
library(tidyverse)
library(R.utils)

setwd("C:/Users/tumok/Documents/Projects/AnalyseSpatialTranscriptomics")

getwd()
list.files("./processed")


gzip ("./processed/matrix.mtx")
gzip ("./processed/barcodes.tsv")
gzip ("./processed/features.tsv")

Raw_data <- Read10X(data.dir = "./processed", gene.column = 1)
metadata <- read.csv("./processed/metadata.csv")

rownames(metadata) <- metadata[["X"]]
metadata[["X"]]<-NULL

processedSeuratObject <- CreateSeuratObject(counts = Raw_data, meta.data = metadata)

view(processedSeuratObject@meta.data)

saveRDS(processedSeuratObject, file="./processed/processedSeuratObject.RDS")

######################################################################################
processedSeuratObject <- readRDS("./processed/processedSeuratObject.RDS")
processedSeuratObject.list <- SplitObject(processedSeuratObject, split.by = "cell_type")
