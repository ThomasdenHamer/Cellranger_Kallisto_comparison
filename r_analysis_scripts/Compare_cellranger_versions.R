library(tidyverse)
library(Seurat)
library(hdf5r)
library(patchwork)

CELLRANGERFOL="/exports/sascstudent/thomas/output/cellranger/"
#list.files(paste0(CELLRANGERFOL, "3kpbmc_1.1.0/filtered_gene_bc_matrices/hg19/"))

{
  cell_1.data <- Read10X(paste0(CELLRANGERFOL, "3kpbmc_1.1.0/filtered_gene_bc_matrices/hg19/"))
  cell_1 <- CreateSeuratObject(counts = cell_1.data, project = "1.1.0")
  remove(cell_1.data)
  
  cell_7_intron.data <- Read10X(paste0(CELLRANGERFOL, "3kpbmc_7.0.1_intron/3kpbmc/outs/filtered_feature_bc_matrix/"))
  cell_7_intron <- CreateSeuratObject(counts = cell_7_intron.data, project = "7.0.1_intron")
  remove(cell_7_intron.data)
  
  cell_7_nointron.data <- Read10X(paste0(CELLRANGERFOL, "3kpbmc_7.0.1_nointron/3kpbmc/outs/filtered_feature_bc_matrix/"))
  cell_7_nointron <- CreateSeuratObject(counts = cell_7_nointron.data, project = "7.0.1_nointron")
  remove(cell_7_nointron.data)
  
  cell_6_intron.data <- Read10X(paste0(CELLRANGERFOL, "3kpbmc-intron-6.1.2/3kpbmc/outs/filtered_feature_bc_matrix/"))
  cell_6_intron <- CreateSeuratObject(counts = cell_6_intron.data, project = "6.1.2_intron")
  remove(cell_6_intron.data)
  
  cell_6_nointron.data <- Read10X(paste0(CELLRANGERFOL, "3kpbmc-nointron-6.1.2/3kpbmc/outs/filtered_feature_bc_matrix/"))
  cell_6_nointron <- CreateSeuratObject(counts = cell_6_nointron.data, project = "6.1.2_nointron")
  remove(cell_6_nointron.data)
  gc()
}

#ggsave("filename.png", plot = x, path = "exports/sascstudent/thomas/R_output/cellranger_version_compare/", device = "png")

cell_1
cell_7_nointron
cell_7_intron

cell_1[["percent.mt"]] <- PercentageFeatureSet(cell_1, pattern = "^MT-")
cell_7_nointron[["percent.mt"]] <- PercentageFeatureSet(cell_7_nointron, pattern = "^MT-")
cell_7_intron[["percent.mt"]] <- PercentageFeatureSet(cell_7_intron, pattern = "^MT-")
cell_6_nointron[["percent.mt"]] <- PercentageFeatureSet(cell_6_nointron, pattern = "^MT-")
cell_6_intron[["percent.mt"]] <- PercentageFeatureSet(cell_6_intron, pattern = "^MT-")


cell_1.vln <- VlnPlot(cell_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")) + ylim(0,20)
cell_7_nointron.vln <- VlnPlot(cell_7_nointron, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")) + ylim(0,20)
cell_7_intron.vln <- VlnPlot(cell_7_intron, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")) + ylim(0,20)
VlnPlot(cell_6_nointron, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")) + ylim(0,20)
VlnPlot(cell_6_intron, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")) + ylim(0,20)

cell_1.vln
cell_7_nointron.vln
cell_7_intron.vln

cell_1.vln / cell_7_nointron.vln / cell_7_intron.vln
remove(cell_1.vln, cell_7_intron.vln, cell_7_nointron.vln)
gc()

cell_1 <- subset(cell_1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
cell_7_nointron <- subset(cell_7_nointron, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
cell_7_intron <- subset(cell_7_intron, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
cell_6_nointron <- subset(cell_6_nointron, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
cell_6_intron <- subset(cell_6_intron, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cell_1
cell_7_nointron
cell_7_intron
cell_6_nointron
cell_6_intron

cell_1 <- NormalizeData(cell_1)
cell_7_nointron <- NormalizeData(cell_7_nointron)
cell_7_intron <- NormalizeData(cell_7_intron)
cell_6_nointron <- NormalizeData(cell_6_nointron)
cell_6_intron <- NormalizeData(cell_6_intron)
###########

cell_1 <- FindVariableFeatures(cell_1, selection.method = "vst", nfeatures = 2000)
cell_7_nointron <- FindVariableFeatures(cell_7_nointron, selection.method = "vst", nfeatures = 2000)
cell_7_intron <- FindVariableFeatures(cell_7_intron, selection.method = "vst", nfeatures = 2000)
cell_6_nointron <- FindVariableFeatures(cell_6_nointron, selection.method = "vst", nfeatures = 2000)
cell_6_intron <- FindVariableFeatures(cell_6_intron, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(cell_1), 10)
cell_1_varfeat <- VariableFeaturePlot(cell_1)
cell_1_varfeat <- LabelPoints(plot = cell_1_varfeat, points = top10, repel = TRUE)
cell_1_varfeat
top10

top10 <- head(VariableFeatures(cell_7_nointron), 10)
cell_7_nointron_varfeat <- cell_7_nointron %>% VariableFeaturePlot(.)
cell_7_nointron_varfeat <- LabelPoints(cell_7_nointron_varfeat, points = top10, repel = TRUE)
cell_7_nointron_varfeat
top10

top10 <- head(VariableFeatures(cell_7_intron), 10)
cell_7_intron_varfeat <- VariableFeaturePlot(cell_7_intron)
cell_7_intron_varfeat <- LabelPoints(cell_7_intron_varfeat, points = top10, repel = T)
cell_7_intron_varfeat
top10

remove(cell_1_varfeat, cell_7_intron_varfeat, cell_7_nointron_varfeat, top10)
gc()

cell_1.genes <- rownames(cell_1)
cell_1 <- ScaleData(cell_1, features = cell_1.genes)

cell_7_nointron.genes <- rownames(cell_7_nointron)
cell_7_nointron <- ScaleData(cell_7_nointron, features = cell_7_nointron.genes)

cell_7_intron.genes <- rownames(cell_7_intron)
cell_7_intron <- ScaleData(cell_7_intron, features = cell_7_intron.genes)

cell_6_nointron.genes <- rownames(cell_6_nointron)
cell_6_nointron <- ScaleData(cell_6_nointron, features = cell_6_nointron.genes)

cell_6_intron.genes <- rownames(cell_6_intron)
cell_6_intron <- ScaleData(cell_6_intron, features = cell_6_intron.genes)


cell_1 <- RunPCA(cell_1, features = VariableFeatures(cell_1))
cell_7_nointron <- RunPCA(cell_7_nointron, features = VariableFeatures(cell_7_nointron))
cell_7_intron <- RunPCA(cell_7_intron, features = VariableFeatures(cell_7_intron))
cell_6_nointron <- RunPCA(cell_6_nointron, features = VariableFeatures(cell_6_nointron))
cell_6_intron <- RunPCA(cell_6_intron, features = VariableFeatures(cell_6_intron))


ElbowPlot(cell_1)
ElbowPlot(cell_7_nointron)
ElbowPlot(cell_7_intron)
ElbowPlot(cell_6_nointron)
ElbowPlot(cell_6_intron)

cell_1 <- FindNeighbors(cell_1, dims = 1:10)
cell_1 <- FindClusters(cell_1, resolution = 0.5)

cell_7_nointron <- FindNeighbors(cell_7_nointron, dims = 1:10)
cell_7_nointron <- FindClusters(cell_7_nointron, resolution = 0.5)

cell_7_intron <- FindNeighbors(cell_7_intron, dims = 1:10)
cell_7_intron <- FindClusters(cell_7_intron, resolution = 0.5)

cell_6_nointron <- FindNeighbors(cell_6_nointron, dims=1:10)
cell_6_nointron <- FindClusters(cell_6_nointron, resolution = 0.5)

cell_6_intron <- FindNeighbors(cell_6_intron, dims=1:10)
cell_6_intron <- FindClusters(cell_6_intron, resolution = 0.5)


cell_1 <- RunUMAP(cell_1, dims = 1:10)
cell_7_nointron <- RunUMAP(cell_7_nointron, dims = 1:10)
cell_7_intron <- RunUMAP(cell_7_intron, dims=1:10)
cell_6_nointron <- RunUMAP(cell_6_nointron, dims=1:10)
cell_6_intron <- RunUMAP(cell_6_intron, dims=1:10)

DimPlot(cell_1, reduction = "umap")
DimPlot(cell_7_nointron, reduction = 'umap')
DimPlot(cell_7_intron, reduction = 'umap')
DimPlot(cell_6_nointron, reduction = 'umap')
DimPlot(cell_6_intron, reduction = 'umap')

cell_1.varfeat <- head(VariableFeatures(cell_1), 9)
FeaturePlot(cell_1, features = cell_1.varfeat)
FeaturePlot(cell_7_nointron, features = cell_1.varfeat)
FeaturePlot(cell_7_intron, features = cell_1.varfeat)
FeaturePlot(cell_6_nointron, features = cell_1.varfeat)
FeaturePlot(cell_6_intron, features = cell_1.varfeat)

#IGLL5


# get_cluster_count_backup <- function(cell.cluster.num, cell.cluster.name, other.projectname) {
#   print(cell.cluster.name)
#   cell.cluster.barcodes <- cell_7_nointron[[]] %>%
#     filter(seurat_clusters == cell.cluster.num) %>%
#     rownames()
#   test <- cell_7_intron[[]] %>%
#     select(seurat_clusters) %>%
#     rownames_to_column() %>%
#     filter(rowname %in% cell.cluster.barcodes) %>%
#     column_to_rownames() %>%
#     group_by(seurat_clusters) %>%
#     summarise(cell.cluster.name = length(seurat_clusters)) %>%
#     column_to_rownames("seurat_clusters")
#   rownames(test) <- lapply(rownames(test), function(x) paste0(other.projectname, ".clust_", x))
#   colnames(test) <- cell.cluster.name
#   test <- test %>% rownames_to_column()
#   return(test)
# }

#cell.cluster.num = number of the ckuster (0 through 8)
#cell.cluster.name = 7.0.1_nointron.clust_0 = projectname + clust_ + cell.cluster.num

#Clusters similarity matrix of cell counts
get_cluster_count <- function(cell.cluster.num, project1, project2) {
  cell.cluster.barcodes <- project1[[]] %>%
    filter(seurat_clusters == cell.cluster.num) %>%
    rownames()
  project1.cluster.name <- paste0(Project(project1),".clust_",cell.cluster.num)
  test <- project2[[]] %>%
    select(seurat_clusters) %>%
    rownames_to_column() %>%
    filter(rowname %in% cell.cluster.barcodes) %>%
    column_to_rownames() %>%
    group_by(seurat_clusters) %>%
    summarise(project1.cluster.name = length(seurat_clusters)) %>%
    column_to_rownames("seurat_clusters")
  rownames(test) <- lapply(rownames(test), function(x) paste0(Project(project2), ".clust_", x))
  colnames(test) <- paste0(Project(project1),".clust_",cell.cluster.num)
  test <- test %>% rownames_to_column()
  return(test)
}

get_shared_cells <- function(project1, project2){
  clus_list <- lapply(levels(project1@meta.data$seurat_clusters), FUN = function(x) get_cluster_count(x, project1, project2))
  clust_share <- data.frame(x=0:8)
  rownames(clust_share) <- lapply(seq(1:9)-1, function(x) paste0(Project(project2), ".clust_", x))
  clust_share <- clust_share %>% rownames_to_column()
  clust_share <- subset(clust_share, select = -c(x))
  
  for (i in seq(1:length(clus_list))){
    clust_share <- left_join(clust_share, clus_list[[i]], by="rowname")
  }
  clust_share <- clust_share %>%
    column_to_rownames() %>%
    mutate_all(~replace(., is.na(.), 0))
  for (i in seq(1:ncol(clust_share))){
    names(clust_share)[names(clust_share) == colnames(clust_share)[i]] <-paste0(colnames(clust_share)[i], ".", as.character(sum(clust_share[i])))
  }
  for (i in seq(1:nrow(clust_share))){
    rownames(clust_share)[rownames(clust_share) == rownames(clust_share)[i]] <- paste0(rownames(clust_share)[i], ".", sum(clust_share[i,]))
  }
  remove(i, clus_list)
  return(clust_share)
}

#clust_share <- get_shared_cells(cell_1, cell_6_nointron)
#clust_share <- get_shared_cells(cell_1, cell_7_nointron)
#clust_share <- get_shared_cells(cell_1, cell_6_intron)
#clust_share <- get_shared_cells(cell_1, cell_7_intron)

clust_share <- get_shared_cells(cell_7_intron, cell_6_intron)
clust_share <- get_shared_cells(cell_7_nointron, cell_6_nointron)
clust_share <- get_shared_cells(cell_7_intron, cell_6_nointron)
clust_share <- get_shared_cells(cell_7_nointron, cell_6_intron)
#clust_share

################################################################################
#Clusters similarity matrix percentage
remove(clust_share.perc)
{
  clust_share.perc <- clust_share
  for (i in seq(1:length(colnames(clust_share)))) {
    newcol <- clust_share[,i] / sum(clust_share[,i])
    clust_share.perc[,i] <- round(newcol * 100, 2)
  }
  remove(i, newcol)
}
#clust_share.perc
#Clusters similarity heatmap
{
  sample1.clusters <- colnames(clust_share)
  sample2.clusters <- rownames(clust_share)
  data <- expand_grid(sample1.clusters, sample2.clusters)
  for (i in seq(1:nrow(data))){
    data[i, "cells"] <- clust_share.perc[as.character(data[i,"sample2.clusters"]), as.character(data[i,"sample1.clusters"])]
  }
  remove(sample1.clusters,sample2.clusters,i)
}
ggplot(data, aes(sample1.clusters, sample2.clusters, fill = cells)) + 
  geom_tile() + 
  geom_text(aes(label = cells), color = "#FFFFFF") +
  theme(axis.text.x = element_text(angle=70, vjust = 0.5)) + 
  scale_fill_gradient(breaks=c(100,75,50,25,0), labels=c("100%", "75%", "50%", "25%","0%"))


cell_1

cell_6_intron
cell_7_intron

cell_6_nointron
cell_7_nointron



