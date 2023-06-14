library(tidyverse)
library(Seurat)
library(future)
library(SeuratDisk)
library(batchelor)
library(SeuratWrappers)
library(gprofiler2)

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

kal_inhouse
kal_inhouse <- NormalizeData(kal_inhouse, normalization.method = "LogNormalize", scale.factor = 100000)
kal_inhouse <- FindVariableFeatures(kal_inhouse)
kal_inhouse <- ScaleData(kal_inhouse, features = rownames(kal_inhouse))
kal_inhouse <- RunPCA(kal_inhouse)
Idents(kal_inhouse)
DimPlot(kal_inhouse, reduction = "pca", group.by = "sample_name")
DimPlot(kal_inhouse, reduction = "pca", group.by = "sample_name", 
        cells = WhichCells(kal_inhouse, idents = c("0","1")))





DimPlot(cell_inhouse, reduction = "umap", group.by = "seurat_clusters")
DimPlot(kal_inhouse, reduction = "umap", group.by = "seurat_clusters")

cell_inhouse[["combined_orig_clusters"]] <- cell_inhouse[["seurat_clusters"]]
kal_inhouse[["combined_orig_clusters"]] <- kal_inhouse[["seurat_clusters"]]

cell_inhouse[["cell_orig_clusters"]] <- cell_inhouse[["seurat_clusters"]]
kal_inhouse[["kal_orig_clusters"]] <- kal_inhouse[["seurat_clusters"]]

Idents(cell_inhouse) <- "seurat_clusters"
DimPlot(cell_inhouse, reduction = "umap", group.by = "seurat_clusters", cells = WhichCells(cell_inhouse, idents = c("0","1")))

Idents(kal_inhouse) <- "seurat_clusters"
DimPlot(kal_inhouse, reduction = "umap", group.by = "seurat_clusters", cells = WhichCells(kal_inhouse, idents = c("0","1")))

kal_inhouse[["origin"]] <- "kallisto"
cell_inhouse[["origin"]] <- "cellranger"


{
Merged_complete <- merge(x = cell_inhouse, y = kal_inhouse, add.cell.ids = c("cell", "kal"))
Merged_complete <- RunFastMNN(object.list = SplitObject(Merged_complete, split.by = "sample_name"))
DimPlot(Merged_complete, reduction = "mnn", group.by = "sample_name")
Merged_complete <- FindNeighbors(Merged_complete, c(1:15), reduction = "mnn")
Merged_complete <- FindClusters(Merged_complete, resolution = 0.27)
Merged_complete <- RunUMAP(Merged_complete, reduction = "mnn", dims = c(1:15))
}
DimPlot(Merged_complete, reduction = "umap", group.by = "seurat_clusters")
DimPlot(Merged_complete, reduction = "umap", group.by = "sample_name")
DimPlot(Merged_complete, reduction = "umap", group.by = "origin")
DimPlot(Merged_complete, reduction = "umap", group.by = "combined_orig_clusters")

saveRDS(Merged_complete, file = "/exports/sascstudent/thomas/R_output/merged_inhouse.rds")

#Save Merged to H5ad object
Merged_complete[["merged_clusters_ch"]] <- clust_char(Merged_complete)
SaveH5Seurat(Merged_complete, filename = "/exports/sascstudent/thomas/R_output/merged_complete.h5seurat", overwrite = T)
Convert("/exports/sascstudent/thomas/R_output/merged_complete.h5seurat", dest = "h5ad", overwrite = T)



#Subsetting the datasets to only include cells from cluster 0 and 1
{
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
}
cell_inhouse_cl01 <- NormalizeData(cell_inhouse_cl01, normalization.method = "LogNormalize", scale.factor = 100000)
cell_inhouse_cl01 <- FindVariableFeatures(cell_inhouse_cl01)
cell_inhouse_cl01 <- ScaleData(cell_inhouse_cl01, features = rownames(cell_inhouse_cl01))
cell_inhouse_cl01 <- RunPCA(cell_inhouse_cl01, features = VariableFeatures(object = cell_inhouse_cl01))
cell_inhouse_cl01 <- CellCycleScoring(cell_inhouse_cl01, s.features = s.genes, g2m.features = g2m.genes, set.indent=TRUE)


#Kallisto
{
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
kal_inhouse_cl01[["origin"]] <- "Kallisto"
kal_inhouse_cl01[["seurat_clusters_kal_orig"]] <- kal_inhouse_cl01$seurat_orig_clusters
gc()
}
kal_inhouse_cl01 <- NormalizeData(kal_inhouse_cl01, normalization.method = "LogNormalize", scale.factor = 100000)
kal_inhouse_cl01 <- FindVariableFeatures(kal_inhouse_cl01)
kal_inhouse_cl01 <- ScaleData(kal_inhouse_cl01, features = rownames(kal_inhouse_cl01))
kal_inhouse_cl01 <- RunPCA(kal_inhouse_cl01, features = VariableFeatures(object = kal_inhouse_cl01))

table(kal_inhouse_cl01[["sample_name"]])


#Merged

{
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
}

Merged_cl01 <- NormalizeData(Merged_cl01, normalization.method = "LogNormalize", scale.factor = 100000)
Merged_cl01 <- FindVariableFeatures(Merged_cl01)
Merged_cl01 <- ScaleData(Merged_cl01, features = rownames(Merged_cl01))
Merged_cl01 <- RunPCA(Merged_cl01, features = VariableFeatures(object = Merged_cl01))

table(Merged_cl01[["sample_name"]])

#Markers
Idents(Merged_cl01) <- "origin"
origin_markers <- FindAllMarkers(Merged_cl01, min.pct = 0.25, only.pos = T, 
                                 logfc.threshold = 0.25)
origin_markers <- origin_markers %>% filter(pct.1>0.6 & p_val_adj <0.05 )
origin_markers %>% group_by(cluster) %>% slice_max(n=10, order_by = avg_log2FC)
origin_markers$gene %>% cat(sep = ",")

DoHeatmap(Merged_cl01, features = origin_markers$gene, group.by = "origin")

origin_marker_pc <- as.data.frame(Loadings(Merged_cl01, reduction = "pca")[row.names(Loadings(Merged_cl01, reduction = "pca")) 
                                         %in% origin_markers$gene,1:15])
origin_marker_pc %>% slice_max(n=length(origin_marker_pc$PC_1), order_by = PC_1)

origin_markers %>% filter(pct.2 < 0.2) %>% .$gene
origin_marker_pc[]


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

Idents(Merged_cl01) <- "origin"
cell_kal_markers <- FindAllMarkers(Merged_cl01, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
cell_kal_markers <- cell_kal_markers %>% filter(pct.1>0.6 & p_val_adj <0.05)
markers <- cell_kal_markers %>% filter(pct.1 > 0.25) %>% filter(pct.2 > 0.25) %>% 
  group_by(cluster) %>% slice_max(n=4, order_by = avg_log2FC)
markers

RidgePlot(Merged_cl01, features = markers$gene, ncol = 2, group.by = "origin")
round()
plot <- DotPlot(Merged_cl01, features = markers$gene) + aes(label = round(pct.exp, digits = 2)) + 
  RotatedAxis() + geom_text(nudge_y = -0.25)
plot
DotPlot(Merged_cl01, features = markers$gene, group.by = "sample_name", split.by = "origin") + RotatedAxis()

plot <- DimPlot(Merged_cl01, reduction = "mnn", group.by = "sample_name")
plot
LabelPoints(plot = plot, points = TopCells(object = Merged_cl01[["mnn"]]))



Idents(kal_inhouse_cl01) <- "sample_name"
sample_markers <- FindAllMarkers(kal_inhouse_cl01, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
sample_markers <- FindMarkers(kal_inhouse_cl01, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, 
                              group.by = "sample_name")
sample_markers <- sample_markers %>% filter(pct.1>0.6 & p_val_adj <0.05) %>% 
  filter(pct.1 > 0.25) %>% filter(pct.2 > 0.25)
markers <- sample_markers %>% group_by(cluster) %>% slice_max(n=2, order_by = avg_log2FC)
markers
 
RidgePlot(kal_inhouse_cl01, features = markers$gene, ncol = 2)
DotPlot(kal_inhouse_cl01, features = markers$gene, ncol = 2, group.by = "sample_name")


DimPlot(kal_inhouse_cl01, reduction = "mnn", group.by = "sample_name")

top_mnn_merged <- TopCells(kal_inhouse_cl01[["mnn"]], ncells = 400)
DimPlot(kal_inhouse_cl01, reduction = "mnn", cells.highlight = top_mnn_merged)

RidgePlot(kal_inhouse_cl01, features = markers$gene, ncol = 2)

DimPlot(cell_inhouse_cl01, reduction = "umap", group.by = "sample_name")
DimPlot(kal_inhouse_cl01, reduction = "umap", group.by = "sample_name")


#Complete merged dataset markers
DimPlot(Merged_complete, reduction = "umap", group.by = "seurat_clusters")
Idents(Merged_complete) <- "seurat_clusters"

#Cellranger POSITIVE markers
Merged_markers_3_4 <- FindMarkers(Merged_complete, ident.1 = "3", ident.2 = "4",
                                       min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
Merged_markers_3_4 %>% filter(pct.1 > 0.6, p_val_adj < 0.05)
Merged_markers_3_4 %>% filter(pct.1 > 0.6, p_val_adj < 0.05) %>% rownames() %>% cat(sep=",")

gostres <- gost(query = Merged_markers_3_4 %>% filter(pct.1 > 0.6, p_val_adj < 0.05) %>% rownames(), 
                organism = "hsapiens", user_threshold = 0.5)
gostplot(gostres)
length(rownames(gostres$result))

#Kallisto POSITIVE markers
Merged_markers_4_3 <- FindMarkers(Merged_complete, ident.1 = "4", ident.2 = "3",
                                  min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
Merged_markers_4_3 %>% filter(pct.1 > 0.6, p_val_adj < 0.05)
Merged_markers_4_3 %>% filter(pct.1 > 0.6, p_val_adj < 0.05) %>% rownames() %>% cat(sep=",")

gostres <- gost(query = Merged_markers_4_3 %>% filter(pct.1 > 0.6, p_val_adj < 0.05) %>% rownames(), 
                organism = "hsapiens", user_threshold = 0.5)
gostplot(gostres)
length(rownames(gostres$result))

####################
#Cellranger POSITIVE markers
Merged_markers_3_4 <- FindMarkers(Merged_complete, ident.1 = "3", ident.2 = "4",
                                  min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
Merged_markers_3_4 %>% filter(pct.1 > 0.6, p_val_adj < 0.05)
Merged_markers_3_4 %>% filter(pct.1 > 0.6, p_val_adj < 0.05) %>% rownames() %>% cat(sep=",")

gostres <- gost(query = Merged_markers_3_4 %>% filter(pct.1 > 0.6, p_val_adj < 0.05) %>% rownames(), 
                organism = "hsapiens", user_threshold = 0.5)
gostplot(gostres)
length(rownames(gostres$result))

#Kallisto POSITIVE markers
Merged_markers_4_3 <- FindMarkers(Merged_complete, ident.1 = "4", ident.2 = "3",
                                  min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
Merged_markers_4_3 %>% filter(pct.1 > 0.6, p_val_adj < 0.05)
Merged_markers_4_3 %>% filter(pct.1 > 0.6, p_val_adj < 0.05) %>% rownames() %>% cat(sep=",")

gostres <- gost(query = Merged_markers_4_3 %>% filter(pct.1 > 0.6, p_val_adj < 0.05) %>% rownames(), 
                organism = "hsapiens", user_threshold = 0.5)
gostplot(gostres)
length(rownames(gostres$result))

###########################

hPSC_merged <- readRDS("/exports/sascstudent/thomas/R_output/merged_inhouse.rds")
hPSC_merged
View(hPSC_merged@meta.data)
DimPlot(hPSC_merged, reduction = "umap", group.by = "origin")
DimPlot(hPSC_merged, reduction = "umap", group.by = "cell_orig_clusters")



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




