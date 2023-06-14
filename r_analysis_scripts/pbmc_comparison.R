library(tidyverse)
library(Seurat)
library(Matrix)
library(batchelor)
library(SeuratWrappers)
library(gprofiler2)

library(future)
library(SeuratDisk)



CELLRANGERFOL="/exports/sascstudent/thomas/output/cellranger/"
KALBUSFOL="/exports/sascstudent/thomas/output/bustools/"

clust_char <- function(seuratobj, metafeature){
  outputdf <- seuratobj[[metafeature]] %>% rownames_to_column()
  factor_columns <- sapply(outputdf, is.factor)
  outputdf[factor_columns] <- lapply(outputdf[factor_columns], as.character)
  outputdf <- outputdf %>% column_to_rownames()
  return(outputdf)
}

# Slightly modified from BUSpaRse, just to avoid installing a few dependencies not used here
read_count_output <- function(dir, name) {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "dgCMatrix")
  # The matrix read has cells in rows
  ge <- ".genes.txt"
  genes <- readLines(file(paste0(dir, "/", name, ge)))
  barcodes <- readLines(file(paste0(dir, "/", name, ".barcodes.txt")))
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}

tr2g <- read_tsv(paste0("/exports/sascstudent/thomas/ref_indices/GRch38_kallisto_cdna", "/tr2g.txt"), progress = F,
                 col_names = c("transcript", "gene", "gene_symbol","gene_full_symbol","chr","start","end","strand")) %>% select(-c(transcript, gene_full_symbol, chr, start, end, strand)) %>% distinct()

{
  cell.data <- Read10X("/exports/sascstudent/thomas/output/cellranger/3kpbmc_1.1.0/filtered_gene_bc_matrices/hg19/", strip.suffix = T)
  cell_pbmc_before <- CreateSeuratObject(cell.data, project = "cellranger", min.cells = 3, min.features = 200)
  remove(cell.data)
  
  res_mat <- read_count_output(paste0(KALBUSFOL, "/3k_pbmc_default/counts_unfiltered"), name = "cells_x_genes")
  rownames(res_mat) <- tr2g$gene_symbol[match(rownames(res_mat), tr2g$gene)]
  kal_pbmc_before <- CreateSeuratObject(counts = res_mat, project = "kal_default", min.cells = 3, min.features = 100)
  remove(res_mat)
}


#Kallisto default workflow
kal_pbmc[["percent.mt"]] <- PercentageFeatureSet(kal_pbmc, pattern = "^MT-")
VlnPlot(kal_pbmc, features = c("nCount_RNA","nFeature_RNA","percent.mt"), ncol = 3)
plot1 <- FeatureScatter(kal_pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
plot2 <- FeatureScatter(kal_pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
plot1 | plot2
remove(plot1, plot2)
kal_pbmc <- subset(kal_pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
kal_pbmc
kal_pbmc <- NormalizeData(kal_pbmc)
kal_pbmc <- FindVariableFeatures(kal_pbmc, nfeatures = 2000, selection.method = "vst")
kal_pbmc <- ScaleData(kal_pbmc, features = rownames(kal_pbmc))
kal_pbmc <- RunPCA(kal_pbmc, features = VariableFeatures(kal_pbmc))
ElbowPlot(kal_pbmc)

kal_pbmc <- FindNeighbors(kal_pbmc, dims = 1:10)
kal_pbmc <- FindClusters(kal_pbmc, resolution = 0.5)
kal_pbmc <- RunUMAP(kal_pbmc, dims = 1:10)
DimPlot(kal_pbmc, reduction="umap")


#Cellranger default workflow
cell_pbmc[["percent.mt"]] <- PercentageFeatureSet(cell_pbmc, pattern = "^MT-")
VlnPlot(cell_pbmc, features = c("nCount_RNA","nFeature_RNA","percent.mt"), ncol = 3)
plot1 <- FeatureScatter(cell_pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
plot2 <- FeatureScatter(cell_pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
plot1 | plot2
remove(plot1, plot2)
cell_pbmc <- subset(cell_pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
cell_pbmc
cell_pbmc <- NormalizeData(cell_pbmc)
cell_pbmc <- FindVariableFeatures(cell_pbmc, nfeatures = 2000, selection.method = "vst")
cell_pbmc <- ScaleData(cell_pbmc, features = rownames(cell_pbmc))
cell_pbmc <- RunPCA(cell_pbmc, features = VariableFeatures(cell_pbmc))
ElbowPlot(cell_pbmc)

cell_pbmc <- FindNeighbors(cell_pbmc, dims = 1:10)
cell_pbmc <- FindClusters(cell_pbmc, resolution = 0.5)
cell_pbmc <- RunUMAP(cell_pbmc, dims = 1:10)
DimPlot(cell_pbmc, reduction="umap")

#saveRDS(kal_pbmc, file = "/exports/sascstudent/thomas/R_output/kal_pbmc.rds")
#saveRDS(cell_pbmc, file = "/exports/sascstudent/thomas/R_output/cell_pbmc.rds")

kal_pbmc <- readRDS("/exports/sascstudent/thomas/R_output/kal_pbmc.rds")
cell_pbmc <- readRDS("/exports/sascstudent/thomas/R_output/cell_pbmc.rds")

{
kal_pbmc[["origin"]] <- "kallisto"
cell_pbmc[["origin"]] <- "cellranger"
kal_pbmc[["kallisto_clusters"]] <- clust_char(kal_pbmc)
cell_pbmc[["cellranger_clusters"]] <- clust_char(cell_pbmc)
}


merged_pbmc <- merge(x = cell_pbmc, y = kal_pbmc, add.cell.ids = c("kal","cell"))
merged_pbmc <- subset(merged_pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
merged_pbmc <- NormalizeData(merged_pbmc)
merged_pbmc <- FindVariableFeatures(merged_pbmc, nfeatures = 2000, selection.method = "vst")
merged_pbmc <- ScaleData(merged_pbmc, features = rownames(merged_pbmc))
merged_pbmc <- RunPCA(merged_pbmc, features = VariableFeatures(merged_pbmc))
ElbowPlot(merged_pbmc)

merged_pbmc <- FindNeighbors(merged_pbmc, dims = 1:10)
merged_pbmc <- FindClusters(merged_pbmc, resolution = 0.7)
merged_pbmc <- RunUMAP(merged_pbmc, dims = 1:10)
DimPlot(merged_pbmc, reduction="umap")

#saveRDS(merged_pbmc, file = "/exports/sascstudent/thomas/R_output/merged_pbmc.rds")
merged_pbmc <- readRDS("/exports/sascstudent/thomas/R_output/merged_pbmc.rds")

merged_pbmc[["seurat_clusters_ch"]] <- clust_char(merged_pbmc)
SaveH5Seurat(merged_pbmc, filename = "/exports/sascstudent/thomas/R_output/merged_pbmc.h5seurat", overwrite = T)
Convert("/exports/sascstudent/thomas/R_output/merged_pbmc.h5seurat", dest = "h5ad", overwrite = T)
 

kal_pbmc
cell_pbmc
merged_pbmc

Idents(merged_pbmc) <- "seurat_clusters"

#Kallisto POSITIVE markers
markers_0_1 <- FindMarkers(merged_pbmc, ident.1 = "0", ident.2 = "1", 
                           logfc.threshold = 0.25, min.pct = 0.25, only.pos = T)
markers_0_1 %>% filter(pct.1 > 0.6, p_val_adj < 0.05, pct.2 > 0)
markers_0_1 %>% filter(pct.1 > 0.6, p_val_adj < 0.05, pct.2 > 0) %>% rownames() %>% cat(sep = ",")
markers_0_1 %>% filter(pct.1 > 0.6, p_val_adj < 0.05, pct.2 == 0)
markers_0_1 %>% filter(pct.1 > 0.6, p_val_adj < 0.05, pct.2 == 0) %>% rownames() %>% cat(sep = ",")

markers_0_1.gostres.pct2_pos <- gost(query = markers_0_1 %>% 
                                       filter(pct.1 > 0.6, p_val_adj < 0.05, pct.2 > 0) %>% rownames(),
                                     organism = "hsapiens")
gostplot(markers_0_1.gostres.pct2_pos)

markers_0_1.gostres.pct2_zero <- gost(query = markers_0_1 %>% 
                                       filter(pct.1 > 0.6, p_val_adj < 0.05, pct.2 == 0) %>% rownames(),
                                     organism = "hsapiens")
gostplot(markers_0_1.gostres.pct2_zero)

rownames(markers_0_1) %>% cat(sep = " ")

#Cellranger POSITIVE markers
markers_1_0 <- FindMarkers(merged_pbmc, ident.1 = "1", ident.2 = "0", 
                           logfc.threshold = 0.25, min.pct = 0.25, only.pos = T)
markers_1_0 %>% filter(pct.1 > 0.6, p_val_adj < 0.05, pct.2 > 0)
markers_1_0 %>% filter(pct.1 > 0.6, p_val_adj < 0.05, pct.2 > 0) %>% rownames() %>% cat(sep = ",")
markers_1_0 %>% filter(pct.1 > 0.6, p_val_adj < 0.05, pct.2 == 0)
markers_1_0 %>% filter(pct.1 > 0.6, p_val_adj < 0.05, pct.2 == 0) %>% rownames() %>% cat(sep = ",")


markers_1_0.gostres.pct2_pos <- gost(query = markers_1_0 %>% 
                                       filter(pct.1 > 0.6, p_val_adj < 0.05, pct.2 > 0) %>% rownames(),
                                     organism = "hsapiens")
gostplot(markers_1_0.gostres.pct2_pos)

markers_1_0.gostres.pct2_zero <- gost(query = markers_1_0 %>% 
                                       filter(pct.1 > 0.6, p_val_adj < 0.05, pct.2 == 0) %>% rownames(),
                                     organism = "hsapiens")
gostplot(markers_1_0.gostres.pct2_zero)
rownames(markers_1_0) %>% cat(sep = " ")


#complete dataset Cellranger / Kallisto markers
View(merged_pbmc)
Idents(merged_pbmc) <- "orig.ident"
levels(Idents(merged_pbmc))
#Cellranger POSITIVE markers
cell_markers <- FindMarkers(object = merged_pbmc, ident.1 = "cellranger", ident.2 = "kal_default",
                                logfc.threshold = 0.25, min.pct = 0.25, only.pos = T)

cell_markers %>% filter(pct.1 > 0.6, p_val_adj < 0.05, pct.2 > 0)
cell_markers %>% filter(pct.1 > 0.6, p_val_adj < 0.05, pct.2 > 0) %>% rownames() %>% cat(sep = ",")
cell_markers %>% filter(pct.1 > 0.6, p_val_adj < 0.05, pct.2 == 0)
cell_markers %>% filter(pct.1 > 0.6, p_val_adj < 0.05, pct.2 == 0) %>% rownames() %>% cat(sep = ",")

cell_markers.gostres.pct2_pos <- gost(query = cell_markers %>% 
                                       filter(pct.1 > 0.6, p_val_adj < 0.05, pct.2 > 0) %>% rownames(),
                                     organism = "hsapiens")
gostplot(cell_markers.gostres.pct2_pos)
gost(query = cell_markers %>% filter(pct.1 > 0.6, p_val_adj < 0.05, pct.2 > 0) %>% rownames(),
     organism = "hsapiens", as_short_link = T)
  




#Kallisto POSITIVE markers
kal_markers <- FindMarkers(object = merged_pbmc, ident.1 = "kal_default", ident.2 = "cellranger",
                           logfc.threshold = 0.25, min.pct = 0.25, only.pos = T)





#Saving PBMC to h5ad object
colnames(kal_pbmc[[]])

kal_pbmc[["cellranger_clusters"]] <- clust_char(cell_pbmc, "seurat_clusters")
cell_pbmc[["kallisto_clusters"]] <- clust_char(kal_pbmc, "seurat_clusters")


SaveH5Seurat(kal_pbmc, filename = "/exports/sascstudent/thomas/R_output/kal_pbmc.h5seurat", overwrite = T)
Convert("/exports/sascstudent/thomas/R_output/kal_pbmc.h5seurat", dest = "h5ad", overwrite = T)

SaveH5Seurat(cell_pbmc, filename = "/exports/sascstudent/thomas/R_output/cell_pbmc.h5seurat", overwrite = T)
Convert("/exports/sascstudent/thomas/R_output/cell_pbmc.h5seurat", dest = "h5ad", overwrite = T)







