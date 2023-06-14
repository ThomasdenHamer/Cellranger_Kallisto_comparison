library(tidyverse)
library(Seurat)
library(batchelor)
library(hdf5r)
library(SeuratWrappers)
library(gprofiler2)
library(future)
library(SeuratDisk)

cellranger_path="/exports/sascstudent/thomas/output/cellranger/"


load_cellranger <- function(cellrangerpath, projectname){
  in_counts <- Read10X_h5(cellrangerpath)
  seuratobj <- CreateSeuratObject(counts=in_counts, min.cells=3, min.features=100, project=projectname)
  return(seuratobj)
}


{
  sample_0H_1 <- load_cellranger(paste0(cellranger_path, "sample0H_1/1_Mesonephros/outs/raw_feature_bc_matrix.h5"), "0H_1")
  sample_0H_2 <- load_cellranger(paste0(cellranger_path, "sample0H_2/2_Mesonephros/outs/raw_feature_bc_matrix.h5"), "0H_2")
  sample_48H <- load_cellranger(paste0(cellranger_path, "sample48H/sample48H/outs/raw_feature_bc_matrix.h5"), "48H")
  sample_120H <- load_cellranger(paste0(cellranger_path, "sample120H/sample120H/outs/raw_feature_bc_matrix.h5"), "120H")
  
  cellranger_inhouse_before <- merge(x=sample_0H_1,y=c(sample_0H_2,sample_48H,sample_120H),
                            add.cell.ids=c("batch_2_A","batch_2_B","batch_1_sample_48H","batch_1_sample_120H"))
  remove(sample_0H_1,sample_0H_2,sample_48H,sample_120H)
  
  sample_mapping <- read_tsv("/exports/sascstudent/thomas/hiPSC_susana/sample_name_full_dataset.tsv", 
                             skip = 1, col_names = c("rownames","sample_name")) %>% column_to_rownames(., var = "rownames")
  cellranger_inhouse_before[["sample_name"]] <- sample_mapping
  remove(sample_mapping, cellranger_path)
}

#Parameters
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes





{
  cellranger_inhouse
  cellranger_inhouse[["percent.mt"]] <- PercentageFeatureSet(cellranger_inhouse, pattern = "^MT-")
  cellranger_inhouse <- subset(cellranger_inhouse, subset = nFeature_RNA > 2000 & nFeature_RNA < 7000 & percent.mt < 10 & percent.mt > 0.1 & nCount_RNA < 100000)
  gc(verbose = F)
  cellranger_inhouse <- NormalizeData(cellranger_inhouse, normalization.method = "LogNormalize", scale.factor = 100000)
  cellranger_inhouse
  
  cellranger_inhouse <- FindVariableFeatures(cellranger_inhouse, selection.method = "vst", nfeatures = 2000)
  cellranger_inhouse <- ScaleData(cellranger_inhouse, features = rownames(cellranger_inhouse))
  gc(verbose = F)
  cellranger_inhouse <- CellCycleScoring(cellranger_inhouse, s.features = s.genes, g2m.features = g2m.genes, set.indent=TRUE)
  cellranger_inhouse <- RunPCA(cellranger_inhouse, features = VariableFeatures(object = cellranger_inhouse))
  
  cellranger_inhouse <- RunFastMNN(object.list = SplitObject(cellranger_inhouse, split.by = "sample_name"))
  cellranger_inhouse <- FindNeighbors(cellranger_inhouse, dims = c(1:15), reduction = "mnn")
  cellranger_inhouse <- FindClusters(cellranger_inhouse, resolution = 0.2)
  cellranger_inhouse <- RunUMAP(cellranger_inhouse, reduction = "mnn", dims = c(1:15))
}


saveRDS(cellranger_inhouse, file = "/exports/sascstudent/thomas/R_output/cellranger_inhouse.rds")
DimPlot(cellranger_inhouse, reduction = "umap", group.by = "seurat_clusters")
DimPlot(cellranger_inhouse, reduction = "umap", group.by = "sample_name")


paper_inhouse <- readRDS("/exports/sascstudent/thomas/hiPSC_susana/hiPSC_complete_mnn.rds")
DimPlot(paper_inhouse, group.by = "seurat_clusters")
DimPlot(paper_inhouse, group.by = "sample_name")

rownames(cellranger_inhouse[[]])
rownames(paper_inhouse[[]])
length(rownames(cellranger_inhouse[[]]))
length(rownames(paper_inhouse[[]]))
length(intersect(rownames(cellranger_inhouse[[]]), rownames(paper_inhouse[[]])))

paper_inhouse <- RunFastMNN(object.list = SplitObject(paper_inhouse, split.by = "sample_name"))
paper_inhouse <- FindNeighbors(paper_inhouse, dims = c(1:15), reduction = "mnn")
paper_inhouse <- FindClusters(paper_inhouse, resolution = 0.3)
paper_inhouse <- RunUMAP(paper_inhouse, reduction = "mnn", dims = c(1:15))
DimPlot(paper_inhouse, reduction = "umap", group.by = "seurat_clusters")

#
paper_inhouse
cellranger_inhouse











