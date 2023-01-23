library(tidyverse)
library(Seurat)
library(future)
library(SeuratDisk)
library(batchelor)
library(SeuratWrappers)

clust_char <- function(seuratobj){
  outputdf <- seuratobj[["seurat_clusters"]] %>% rownames_to_column()
  factor_columns <- sapply(outputdf, is.factor)
  outputdf[factor_columns] <- lapply(outputdf[factor_columns], as.character)
  outputdf <- outputdf %>% column_to_rownames()
  return(outputdf)
}

plan("multisession", workers = 4)
plan("sequential")
plan()

cell_inhouse <- readRDS("/exports/sascstudent/thomas/hiPSC_susana/hiPSC_complete_mnn.rds")
kal_inhouse <- readRDS("/exports/sascstudent/thomas/R_output/hiPSC_kallisto_inhouse_mnn.rds")

DimPlot(cell_inhouse, reduction = "umap", group.by = "seurat_clusters")
DimPlot(kal_inhouse, reduction = "umap", group.by = "seurat_clusters")

cell_inhouse[["seurat_orig_clusters_ch"]] <- clust_char(cell_inhouse)
kal_inhouse[["seurat_orig_clusters_ch"]] <- clust_char(kal_inhouse)
cell_inhouse[["seurat_orig_clusters"]] <- cell_inhouse[["seurat_clusters"]]
kal_inhouse[["seurat_orig_clusters"]] <- kal_inhouse[["seurat_clusters"]]

Idents(cell_inhouse) <- "seurat_clusters"
DimPlot(cell_inhouse, reduction = "umap", group.by = "seurat_clusters", cells = WhichCells(cell_inhouse, idents = c("0","1")))

Idents(kal_inhouse) <- "seurat_clusters"
DimPlot(kal_inhouse, reduction = "umap", group.by = "seurat_clusters", cells = WhichCells(kal_inhouse, idents = c("0","1")))


#Subsetting the datasets to only include cells from cluster 0 and 1
cell_inhouse_cl01 <- subset(cell_inhouse, idents = c("0","1"))
cell_inhouse_cl01
cell_inhouse_cl01 <- NormalizeData(cell_inhouse_cl01, normalization.method = "LogNormalize", scale.factor = 100000)
cell_inhouse_cl01 <- FindVariableFeatures(cell_inhouse_cl01)
cell_inhouse_cl01 <- ScaleData(cell_inhouse_cl01, features = rownames(cell_inhouse_cl01))
cell_inhouse_cl01 <- RunPCA(cell_inhouse_cl01, features = VariableFeatures(object = cell_inhouse_cl01))

cell_inhouse_cl01 <- RunFastMNN(object.list = SplitObject(cell_inhouse_cl01, split.by = "sample_name"))
cell_inhouse_cl01 <- FindNeighbors(cell_inhouse_cl01, c(1:15), reduction = "mnn")
cell_inhouse_cl01 <- FindClusters(cell_inhouse_cl01, resolution = 0.12)
cell_inhouse_cl01 <- RunUMAP(cell_inhouse_cl01, reduction = "mnn", dims = c(1:15))
DimPlot(cell_inhouse_cl01, group.by = "seurat_clusters")

cell_inhouse_cl01[["origin"]] <- "Cellranger"
cell_inhouse_cl01[["seurat_clusters_cell_orig"]] <- cell_inhouse_cl01$seurat_orig_clusters
gc()

#Kallisto
kal_inhouse_cl01 <- subset(kal_inhouse, idents = c("0","1"))
kal_inhouse_cl01
kal_inhouse_cl01 <- NormalizeData(kal_inhouse_cl01, normalization.method = "LogNormalize", scale.factor = 100000)
kal_inhouse_cl01 <- FindVariableFeatures(kal_inhouse_cl01)
kal_inhouse_cl01 <- ScaleData(kal_inhouse_cl01, features = rownames(kal_inhouse_cl01))
kal_inhouse_cl01 <- RunPCA(kal_inhouse_cl01, features = VariableFeatures(object = kal_inhouse_cl01))

kal_inhouse_cl01 <- RunFastMNN(object.list = SplitObject(kal_inhouse_cl01, split.by = "sample_name"))
kal_inhouse_cl01 <- FindNeighbors(kal_inhouse_cl01, c(1:15), reduction = "mnn")
kal_inhouse_cl01 <- FindClusters(kal_inhouse_cl01, resolution = 0.1)
kal_inhouse_cl01 <- RunUMAP(kal_inhouse_cl01, reduction = "mnn", dims = c(1:15))
DimPlot(kal_inhouse_cl01, group.by = "seurat_clusters")
DimPlot(kal_inhouse_cl01, group.by = "sample_name")

table(kal_inhouse_cl01[["sample_name"]])


kal_inhouse_cl01[["origin"]] <- "Kallisto"
kal_inhouse_cl01[["seurat_clusters_kal_orig"]] <- kal_inhouse_cl01$seurat_orig_clusters


#Merged
Merged_cl01 <- merge(x = cell_inhouse_cl01, y = kal_inhouse_cl01, add.cell.ids = c("Cellranger", "Kallisto"))
Merged_cl01
Merged_cl01[['HTO']] <- NULL
Merged_cl01[['mnn.reconstructed']] <- NULL

Merged_cl01 <- NormalizeData(Merged_cl01, normalization.method = "LogNormalize", scale.factor = 100000)
Merged_cl01 <- FindVariableFeatures(Merged_cl01)
Merged_cl01 <- ScaleData(Merged_cl01, features = rownames(Merged_cl01))
Merged_cl01 <- RunPCA(Merged_cl01, features = VariableFeatures(object = Merged_cl01))

Merged_cl01 <- RunFastMNN(object.list = SplitObject(Merged_cl01, split.by = "sample_name"))
Merged_cl01 <- FindNeighbors(Merged_cl01, c(1:15), reduction = "mnn")
Merged_cl01 <- FindClusters(Merged_cl01, resolution = 0.27)
Merged_cl01 <- RunUMAP(Merged_cl01, reduction = "mnn", dims = c(1:15))
DimPlot(Merged_cl01, reduction = "umap", group.by = "seurat_clusters")
DimPlot(Merged_cl01, reduction = "umap", group.by = "seurat_orig_clusters")
DimPlot(Merged_cl01, reduction = "umap", group.by = "sample_name")
table(Merged_cl01[["sample_name"]])

#Markers
cell_markers <- FindMarkers(Merged_cl01, ident.1 = "Cellranger", ident.2 = "Kallisto",
                                group.by = "origin", min.pct = 0.25, only.pos = T)
cell_markers %>% filter(pct.1>0.6 & p_val_adj <0.05 ) %>% slice_max(n = 10, order_by = avg_log2FC)

kal_markers <- FindMarkers(Merged_cl01, ident.1 = "Kallisto", ident.2 = "Cellranger",
                            group.by = "origin", min.pct = 0.25, only.pos = T)
kal_markers %>% filter(pct.1>0.6 & p_val_adj <0.05 ) %>% slice_max(n = 10, order_by = avg_log2FC)
kal_filt_markers <- kal_markers %>% filter(pct.1>0.6 & p_val_adj <0.05 )

kal_filt_markers %>% slice_max(n = 10, order_by = avg_log2FC) %>% rownames %>% cat(sep = ",")



kal_only_markers <- FindMarkers(Merged_cl01, ident.1 = "3", ident.2 = "0", 
            group.by = "seurat_clusters", only.pos = T)

kal_only_markers %>% filter(pct.1>0.6 & p_val_adj <0.05 ) %>% 
  slice_max(n = 10, order_by = avg_log2FC) %>% 
  rownames() %>%
  cat(sep = ",")

kal_only_markers["BEX2",]
kal_only_markers %>% filter(rownames(.) == "IRX2")

Idents(kal_inhouse_cl01) <- "sample_name"
kal_inhouse_cl01_markers <- FindMarkers(kal_inhouse_cl01, ident.1 = "iPSC_54_3", 
                                        ident.2 = "iPSC_72_1", logfc.threshold = 0.25, min.pct = 0.25)

kal_inhouse_cl01_markers <- FindAllMarkers(kal_inhouse_cl01, logfc.threshold = 0.25, min.pct = 0.25, only.pos = T)


kal_inhouse_cl01_filt_markers <- kal_inhouse_cl01_markers %>% filter(pct.1>0.6 & p_val_adj <0.05 )
kal_inhouse_cl01_filt_markers %>% group_by(cluster) %>% slice_max(n = 3, order_by = avg_log2FC)
kal_inhouse_cl01_filt_markers %>% slice_max(n = 10, order_by = avg_log2FC)  %>% rownames %>% cat(sep = ",")

kal_inhouse_cl01_markers_2 <- FindMarkers(kal_inhouse_cl01, ident.1 = "iPSC_72_1", 
                                        ident.2 = "iPSC_54_3", logfc.threshold = 0.25, min.pct = 0.25)

kal_inhouse_cl01_filt_markers_2 <- kal_inhouse_cl01_markers_2 %>% filter(pct.1>0.6 & p_val_adj <0.05 )
kal_inhouse_cl01_filt_markers_2 %>% slice_max(n = 10, order_by = avg_log2FC)
kal_inhouse_cl01_filt_markers_2 %>% slice_max(n = 10, order_by = avg_log2FC)  %>% rownames %>% cat(sep = ",")





#Save Merged to H5ad object
Merged_cl01[["seurat_clusters_ch"]] <- clust_char(Merged_cl01)
SaveH5Seurat(Merged_cl01, filename = "/exports/sascstudent/thomas/R_output/merged_cl01.h5seurat", overwrite = T)
Convert("/exports/sascstudent/thomas/R_output/merged_cl01.h5seurat", dest = "h5ad", overwrite = T)

#Save Kallisto only cl0 en cl1 to H5ad
kal_inhouse_cl01[["seurat_clusters_ch"]] <- clust_char(kal_inhouse_cl01)
SaveH5Seurat(kal_inhouse_cl01, filename = "/exports/sascstudent/thomas/R_output/kal_inhouse_cl01.h5seurat", overwrite = T)
Convert("/exports/sascstudent/thomas/R_output/kal_inhouse_cl01.h5seurat", dest = "h5ad", overwrite = T)

#Save Cellranger only cl0 en cl1 to H5ad
cell_inhouse_cl01[["seurat_clusters_ch"]] <- clust_char(cell_inhouse_cl01)
SaveH5Seurat(cell_inhouse_cl01, filename = "/exports/sascstudent/thomas/R_output/cell_inhouse_cl01.h5seurat", overwrite = T)
Convert("/exports/sascstudent/thomas/R_output/cell_inhouse_cl01.h5seurat", dest = "h5ad", overwrite = T)




