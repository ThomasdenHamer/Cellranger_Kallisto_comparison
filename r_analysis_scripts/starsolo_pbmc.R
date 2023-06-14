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

pbmc3k <- load_starsolo(paste0(starsolo_path, "pbmc3k/Solo.out/Gene/raw"), "pbmc3k")


#Parameters
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

{
  pbmc3k
  pbmc3k[["percent.mt"]] <- PercentageFeatureSet(pbmc3k, pattern = "^MT-")
  pbmc3k <- subset(pbmc3k, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  gc(verbose = F)
  pbmc3k <- NormalizeData(pbmc3k, normalization.method = "LogNormalize", scale.factor = 100000)
  pbmc3k
  
  pbmc3k <- FindVariableFeatures(pbmc3k)
  pbmc3k <- ScaleData(pbmc3k, features = rownames(pbmc3k))
  gc(verbose = F)
  pbmc3k <- CellCycleScoring(pbmc3k, s.features = s.genes, g2m.features = g2m.genes, set.indent=TRUE)
  pbmc3k <- RunPCA(pbmc3k, features = VariableFeatures(object = pbmc3k))
  
  pbmc3k <- FindNeighbors(pbmc3k, dims = c(1:15))
  pbmc3k <- FindClusters(pbmc3k, resolution = 0.5)
  pbmc3k <- RunUMAP(pbmc3k, dims = c(1:15))
  DimPlot(pbmc3k, reduction = "umap", group.by = "seurat_clusters")
}

saveRDS(pbmc3k, "/exports/sascstudent/thomas/R_output/starsolo_pbmc3k.rds")
#pbmc3k <- readRDS("/exports/sascstudent/thomas/R_output/starsolo_inhouse.rds")

clust_char <- function(seuratobj, metavariable="seurat_clusters"){
  outputdf <- seuratobj[[metavariable]] %>% rownames_to_column()
  factor_columns <- sapply(outputdf, is.factor)
  outputdf[factor_columns] <- lapply(outputdf[factor_columns], as.character)
  outputdf <- outputdf %>% column_to_rownames()
  return(outputdf)
}

pbmc3k[["seurat_clusters_ch"]] <- clust_char(pbmc3k)
SaveH5Seurat(pbmc3k, filename = "/exports/sascstudent/thomas/R_output/starsolo_pbmc3k.h5seurat", overwrite = T)
Convert("/exports/sascstudent/thomas/R_output/starsolo_pbmc3k.h5seurat", dest = "h5ad", overwrite = T)

