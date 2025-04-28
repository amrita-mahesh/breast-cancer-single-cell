library(Seurat)
library(dplyr)
library(magrittr)

#Code to load the dataset as a Seurat object and calculate/plot QC metrics

setwd("D:/Projects/Github side projects/scRNAseq analysis/Data/")

#Load input data as a Seurat object
expression_matrix <- ReadMtx(mtx = "matrix.mtx.gz", cells = "barcodes.tsv.gz", features = "features.tsv.gz", feature.column = 1)
seurat_object <- CreateSeuratObject(expression_matrix, min.cells = 3, min.features = 200)

expression_matrix <- NULL

#Find the percentage of mitochondrial and ribosomal genes in all cells
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
seurat_object[["percent.rb"]] <- PercentageFeatureSet(seurat_object, pattern = "^RP[SL]")

#Extract metadata and add new metrics
metadata <- seurat_object@meta.data
metadata <- metadata %>% rename(nUMI = nCount_RNA, nGene = nFeature_RNA, Sample = orig.ident)
metadata$cells <- rownames(metadata)
metadata$log10GenesPerUMI <- log10(metadata$nGene)/log10(metadata$nUMI)
metadata$project <- rep("BC_project", times = nrow(metadata))

#Map new metadata back to old object
seurat_object@meta.data <- metadata

#Set identities to project name for easier visualization
seurat_object <- SetIdent(seurat_object, value=seurat_object@meta.data$project)

#save(seurat_object, file="D:/Projects/Github side projects/scRNAseq analysis/seurat_object.Rdata")

#Plot the calculated QC metrics
VlnPlot(seurat_object, features = c("nGene", "nUMI"))
VlnPlot(seurat_object, features = c("log10GenesPerUMI"))
VlnPlot(seurat_object, features = c("percent.mt", "percent.rb"))

#Note: Doublet detection and empty drop identification has not been performed on this dataset





