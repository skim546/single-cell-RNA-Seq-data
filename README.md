# single-cell-RNA-Seq-data
Transcriptomic characterisation of haematopoietic stem and progenitor cells from human adult bone marrow, spleen and peripheral blood

library(Seurat)
library(SeuratDisk)
library(Matrix)
library(tidyverse)

input_file <- "~/Desktop/20210204_S_SUBS4_4x_annotated_donors_UMI_counts_HCA.h5ad"
output_file <- "~/Desktop/20210204_S_SUBS4_4x_annotated_donors_UMI_counts_HCA.h5Seurat"
Convert(input_file, dest = "h5seurat", overwrite = TRUE)

seurat_obj <- LoadH5Seurat(output_file)

DefaultAssay(seurat_obj) <- "RNA"

seurat_obj[["nFeature_RNA"]] <- Matrix::colSums(seurat_obj[["RNA"]]@counts > 0)
seurat_obj[["nCount_RNA"]] <- Matrix::colSums(seurat_obj[["RNA"]]@counts)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

head(seurat_obj@meta.data)

VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
