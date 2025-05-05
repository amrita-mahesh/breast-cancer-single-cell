library(Seurat)
library(dplyr)
library(org.Hs.eg.db)


setwd("D:/Projects/Github side projects/scRNAseq analysis/Data/")

#Load input data
expression_matrix <- ReadMtx(mtx = "matrix.mtx.gz", cells = "barcodes.tsv.gz", features = "features.tsv.gz", feature.column = 1)

#Convert ENSEMBL gene IDs to gene symbols
gene_keys <- expression_matrix@Dimnames[[1]]
mapped_genes <- data.frame(mapIds(org.Hs.eg.db, keys = gene_keys, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
mapped_genes$ENSEMBL <- rownames(mapped_genes)
colnames(mapped_genes) <- c("Symbol", "ENSEMBL")
#Accounting for duplicated or unmapped gene symbols
mapped_genes <- mapped_genes %>% 
  mutate(mapped_names = if_else((is.na(Symbol) | duplicated(Symbol)), ENSEMBL, Symbol))

expression_matrix@Dimnames[[1]] <- mapped_genes$mapped_names

#Load Seurat object with new gene names
seurat_object <- CreateSeuratObject(expression_matrix, min.cells = 10, min.features = 200)
expression_matrix <- NULL


seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize")
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)

#Plot the top 20 most variable genes
most_variable <- head(VariableFeatures(seurat_object), 5)
scatter_plot_genes <- VariableFeaturePlot(seurat_object)
LabelPoints(plot = scatter_plot_genes, points = most_variable, repel = TRUE, xnudge = 0, ynudge = 0)


#Scale the data and perform PCA
seurat_object <- ScaleData(object = seurat_object, vars.to.regress = c("nCount_RNA", "orig.ident"))
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
PCAPlot(seurat_object, dims = c(1,2))

#Choose PC cutoff and perform clustering
ElbowPlot(seurat_object)
seurat_object <- FindNeighbors(seurat_object, dims = 1:5)
seurat_object <- FindClusters(seurat_object, resolution = 0.5)

#Run UMAP
seurat_object <- RunUMAP(seurat_object, dims = 1:5, verbose = F)
cells_per_cluster <- table(seurat_object@meta.data$seurat_clusters)
DimPlot(seurat_object,label.size = 4,repel = T,label = T)

