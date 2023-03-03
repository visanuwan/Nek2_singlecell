library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(MAST)
library(dplyr)
library(ggplot2)

#################################################
###### Load and QC data
#################################################

import_h5 <- Read10X_h5(filename = "data/wt_normal/filtered_feature_bc_matrix.h5")
wt_normal <- CreateSeuratObject(counts = import_h5, project = "wt_normal", min.features = 200)
wt_normal[["percent.mt"]] <- PercentageFeatureSet(wt_normal, pattern = "^mt-")
wt_normal <- subset(wt_normal, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)

import_h5 <- Read10X_h5(filename = "data/wt_5T/filtered_feature_bc_matrix.h5")
wt_5T <- CreateSeuratObject(counts = import_h5, project = "wt_5T", min.features = 200)
wt_5T[["percent.mt"]] <- PercentageFeatureSet(wt_5T, pattern = "^mt-")
wt_5T <- subset(wt_5T, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)

import_h5 <- Read10X_h5(filename = "data/ko_normal/filtered_feature_bc_matrix.h5")
ko_normal <- CreateSeuratObject(counts = import_h5, project = "ko_normal", min.features = 200)
ko_normal[["percent.mt"]] <- PercentageFeatureSet(ko_normal, pattern = "^mt-")
ko_normal <- subset(ko_normal, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)

import_h5 <- Read10X_h5(filename = "data/ko_5T/filtered_feature_bc_matrix.h5")
ko_5T <- CreateSeuratObject(counts = import_h5, project = "ko_5T", min.features = 200)
ko_5T[["percent.mt"]] <- PercentageFeatureSet(ko_5T, pattern = "^mt-")
ko_5T <- subset(ko_5T, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)


#################################################
###### Normalize and prepare data for integration
#################################################

nek2.list <- list(wt_normal, wt_5T, ko_normal, ko_5T)
nek2.list <- lapply(X = nek2.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = nek2.list)
nek2.list <- lapply(X = nek2.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#################################################
###### Perform integration
#################################################

nek2.anchors <- FindIntegrationAnchors(object.list = nek2.list, anchor.features = features, reduction = "rpca")
nek2.combined <- IntegrateData(anchorset = nek2.anchors)
DefaultAssay(nek2.combined) <- "integrated"

#################################################
###### Perform the standard workflow for visualization and clustering
#################################################

nek2.combined <- ScaleData(nek2.combined, verbose = FALSE)
nek2.combined <- RunPCA(nek2.combined, npcs = 50, verbose = FALSE)
nek2.combined <- RunUMAP(nek2.combined, reduction = "pca", dims = 1:35)
nek2.combined <- FindNeighbors(nek2.combined, reduction = "pca", dims = 1:35)
nek2.combined <- FindClusters(nek2.combined, resolution = 0.5)

nek2.combined$orig.ident <- factor(x = nek2.combined$orig.ident, levels = c("wt_normal", "ko_normal", "wt_5T", "ko_5T"))

DimPlot(nek2.combined, reduction = "umap",label = T, label.size = 5, repel = F)
DimPlot(nek2.combined, reduction = "umap",label = T, label.size = 5, repel = F, split.by = 'orig.ident') + NoLegend()
DimPlot(nek2.combined, reduction = "umap",label = T, label.size = 5, repel = F, group.by = 'orig.ident') + NoLegend()

#################################################
###### Cell type annotation
#################################################

DefaultAssay(nek2.combined) <- "RNA"

# MM
FeaturePlot(nek2.combined, features = c("Sdc1"))

# Pre-B cells
FeaturePlot(nek2.combined, features = c("Cd74", "H2Ab1"))

# Pro-B cells
FeaturePlot(nek2.combined, features = c("Vpreb3", "Akap12"))

# NK and T cells
FeaturePlot(nek2.combined, features = c("Gzma", "Klra4", "Cd3d", "Cd3e"))

# DCs
FeaturePlot(nek2.combined, features = c("Itgax", "Siglech"))

# Myeloblasts
FeaturePlot(nek2.combined, features = c("Elane", "Mpo", "Ctsg", "Ms4a3"))

# Myelocytes
FeaturePlot(nek2.combined, features = c("Fcnb", "Ltf", "Lcn2"))

# Neutrophils
FeaturePlot(nek2.combined, features = c("Ly6g", "Mmp8", "Cxcr2"))

# Monoblasts
FeaturePlot(nek2.combined, features = c("F13a1", "Irf8", "Ly86"))

# Monocytes
FeaturePlot(nek2.combined, features = c("S100a4", "Pld4", "Csf1r"))

# Macrophages
FeaturePlot(nek2.combined, features = c("Adgre1"))

# Assign cell type to clusters
nek2.combined@active.ident <- factor(nek2.combined@active.ident,
                                     levels = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17'))
nek2.combined$cell_types <- Idents(nek2.combined)
nek2.combined@meta.data <- nek2.combined@meta.data %>% mutate(name_cell_types = cell_types)
levels(nek2.combined$name_cell_types)[levels(nek2.combined$name_cell_types) == '1'] <- 'MM-I'
levels(nek2.combined$name_cell_types)[levels(nek2.combined$name_cell_types) == '2'] <- 'MM-II'
levels(nek2.combined$name_cell_types)[levels(nek2.combined$name_cell_types) == '3'] <- 'MM-III'
levels(nek2.combined$name_cell_types)[levels(nek2.combined$name_cell_types) == '4'] <- 'MM-IV'
levels(nek2.combined$name_cell_types)[levels(nek2.combined$name_cell_types) == '5'] <- 'Pro-B-I'
levels(nek2.combined$name_cell_types)[levels(nek2.combined$name_cell_types) == '6'] <- 'Pro-B-II'
levels(nek2.combined$name_cell_types)[levels(nek2.combined$name_cell_types) == '7'] <- 'Pre-B'
levels(nek2.combined$name_cell_types)[levels(nek2.combined$name_cell_types) == '8'] <- 'NK/T'
levels(nek2.combined$name_cell_types)[levels(nek2.combined$name_cell_types) == '9'] <- 'Neutrophils-I'
levels(nek2.combined$name_cell_types)[levels(nek2.combined$name_cell_types) == '10'] <- 'Neutrophils-II'
levels(nek2.combined$name_cell_types)[levels(nek2.combined$name_cell_types) == '11'] <- 'Myelocytes-I'
levels(nek2.combined$name_cell_types)[levels(nek2.combined$name_cell_types) == '12'] <- 'Myelocytes-II'
levels(nek2.combined$name_cell_types)[levels(nek2.combined$name_cell_types) == '13'] <- 'Myeloblast'
levels(nek2.combined$name_cell_types)[levels(nek2.combined$name_cell_types) == '14'] <- 'Monoblasts'
levels(nek2.combined$name_cell_types)[levels(nek2.combined$name_cell_types) == '15'] <- 'Monocytes'
levels(nek2.combined$name_cell_types)[levels(nek2.combined$name_cell_types) == '16'] <- 'Macrophages'
levels(nek2.combined$name_cell_types)[levels(nek2.combined$name_cell_types) == '17'] <- 'DCs'


#################################################
###### Cell-type specific markers
#################################################

DefaultAssay(nek2.combined) <- "RNA"
all.markers <- FindAllMarkers(object = nek2.combined, test.use = "MAST")
top20.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

#################################################
###### Subcluster analysis
#################################################

# Macrophages
Idents(nek2.combined) <- "name_cell_types"
nek2.combined.mcp <- subset(nek2.combined, idents = "Macrophages")
nek2.combined.mcp <- RunPCA(nek2.combined.mcp, npcs = 50, verbose = TRUE)
nek2.combined.mcp <- FindNeighbors(nek2.combined.mcp, reduction = "pca", dims = 1:50)
nek2.combined.mcp <- FindClusters(nek2.combined.mcp, graph.name = 'integrated_snn',resolution = 1, verbose = TRUE)
nek2.combined.mcp <- RunUMAP(nek2.combined.mcp, reduction = "pca", dims = 1:50, min.dist = 0.3, n.neighbors = 16)
DimPlot(nek2.combined.mcp, reduction = "umap",label = T, label.size = 4, repel = F, split.by = 'orig.ident')

# Macrophages
nek2.combined.nkt <- subset(nek2.combined, idents = "NK/T")
nek2.combined.nkt <- RunPCA(nek2.combined.nkt, npcs = 50, verbose = TRUE)
nek2.combined.nkt <- FindNeighbors(nek2.combined.nkt, reduction = "pca", dims = 1:50)
nek2.combined.nkt <- FindClusters(nek2.combined.nkt, graph.name = 'integrated_snn',resolution = 0.7, verbose = TRUE)
nek2.combined.nkt <- RunUMAP(nek2.combined.nkt, reduction = "pca", dims = 1:50, min.dist = 0.4, n.neighbors = 25)
DimPlot(nek2.combined.nkt, reduction = "umap",label = T, label.size = 4, repel = F, split.by = 'orig.ident')
