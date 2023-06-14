library(tidyverse)
library(batchelor)
library(Seurat)
library(hdf5r)
library(SeuratWrappers)
library(gprofiler2)
library(future)
library(SeuratDisk)

starsolo_path="/exports/sascstudent/thomas/output/STARsolo/"


load_starsolo <- function(starpath, projectname){
  in_counts <- ReadSTARsolo(starpath)
  seuratobj <- CreateSeuratObject(counts=in_counts, min.cells=3, min.features=100, project=projectname)
  return(seuratobj)
}


{
  sample_0H_1 <- load_starsolo(paste0(starsolo_path, "mesonephros_0H_1/Solo.out/Gene/raw"), "0H_1")
  sample_0H_2 <- load_starsolo(paste0(starsolo_path, "mesonephros_0H_2/Solo.out/Gene/raw"), "0H_2")
  sample_48H <- load_starsolo(paste0(starsolo_path, "sample48H/Solo.out/Gene/raw"), "48H")
  sample_120H <- load_starsolo(paste0(starsolo_path, "sample120H/Solo.out/Gene/raw"), "120H")
  
  starsolo_inhouse <- merge(x=sample_0H_1,y=c(sample_0H_2,sample_48H,sample_120H),
        add.cell.ids=c("batch_2_A","batch_2_B","batch_1_sample_48H","batch_1_sample_120H"))
  remove(sample_0H_1,sample_0H_2,sample_48H,sample_120H)
  
  sample_mapping <- read_tsv("/exports/sascstudent/thomas/hiPSC_susana/sample_name_full_dataset.tsv", 
                             skip = 1, col_names = c("rownames","sample_name")) %>% column_to_rownames(., var = "rownames")
  rownames(sample_mapping) <- rownames(sample_mapping) %>% str_sub(end = -3)
  starsolo_inhouse[["sample_name"]] <- sample_mapping
  remove(sample_mapping, starsolo_path)
}

#Parameters
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

{
  starsolo_inhouse
  starsolo_inhouse[["percent.mt"]] <- PercentageFeatureSet(starsolo_inhouse, pattern = "^MT-")
  starsolo_inhouse <- subset(starsolo_inhouse, subset = nFeature_RNA > 2000 & nFeature_RNA < 7000 & percent.mt < 10 & percent.mt > 0.1 & nCount_RNA < 100000)
  gc(verbose = F)
  starsolo_inhouse <- NormalizeData(starsolo_inhouse, normalization.method = "LogNormalize", scale.factor = 100000)
  starsolo_inhouse
  
  starsolo_inhouse <- FindVariableFeatures(starsolo_inhouse)
  starsolo_inhouse <- ScaleData(starsolo_inhouse, features = rownames(starsolo_inhouse))
  gc(verbose = F)
  starsolo_inhouse <- CellCycleScoring(starsolo_inhouse, s.features = s.genes, g2m.features = g2m.genes, set.indent=TRUE)
  starsolo_inhouse <- RunPCA(starsolo_inhouse, features = VariableFeatures(object = starsolo_inhouse))
  
  starsolo_inhouse <- RunFastMNN(object.list = SplitObject(starsolo_inhouse, split.by = "sample_name"))
  starsolo_inhouse <- FindNeighbors(starsolo_inhouse, dims = c(1:15), reduction = "mnn")
  starsolo_inhouse <- FindClusters(starsolo_inhouse, resolution = 0.2)
  starsolo_inhouse <- RunUMAP(starsolo_inhouse, reduction = "mnn", dims = c(1:15))
  DimPlot(starsolo_inhouse, reduction = "umap", group.by = "seurat_clusters")
  DimPlot(starsolo_inhouse, reduction = "umap", group.by = "sample_name")
}

saveRDS(starsolo_inhouse, "/exports/sascstudent/thomas/R_output/starsolo_inhouse.rds")
starsolo_inhouse <- readRDS("/exports/sascstudent/thomas/R_output/starsolo_inhouse.rds")

clust_char <- function(seuratobj, metavariable="seurat_clusters"){
  outputdf <- seuratobj[[metavariable]] %>% rownames_to_column()
  factor_columns <- sapply(outputdf, is.factor)
  outputdf[factor_columns] <- lapply(outputdf[factor_columns], as.character)
  outputdf <- outputdf %>% column_to_rownames()
  return(outputdf)
}

starsolo_inhouse[["seurat_clusters_ch"]] <- clust_char(starsolo_inhouse)
SaveH5Seurat(starsolo_inhouse, filename = "/exports/sascstudent/thomas/R_output/starsolo_inhouse.h5seurat", overwrite = T)
Convert("/exports/sascstudent/thomas/R_output/starsolo_inhouse.h5seurat", dest = "h5ad", overwrite = T)



















