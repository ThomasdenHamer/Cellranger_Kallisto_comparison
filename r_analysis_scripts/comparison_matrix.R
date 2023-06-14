library(tidyverse)
library(Seurat)
library(cluster)
library(pdfCluster, include.only = 'adj.rand.index')
library(fossil, include.only = "rand.index")
library(gprofiler2)
library(patchwork)
library(Matrix)

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
load_starsolo <- function(starpath, projectname){
  in_counts <- ReadSTARsolo(starpath)
  seuratobj <- CreateSeuratObject(counts=in_counts, min.cells=3, min.features=100, project=projectname)
  return(seuratobj)
}

{
  pbmc.data <- Read10X(data.dir = "/exports/sascstudent/thomas/output/cellranger/3kpbmc_7.0.1_nointron/3kpbmc/outs/raw_feature_bc_matrix/")
  cellranger_pbmc_before <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k.cr", min.cells = 3, min.features = 200)
  ####
  pbmc.kalisto.data <-  read_count_output("/exports/sascstudent/thomas/output/bustools/3k_pbmc_default/counts_unfiltered", name = "cells_x_genes")
  tr2g <- read_tsv(paste0("/exports/sascstudent/thomas/ref_indices/GRch38_kallisto_cdna", "/tr2g.txt"), progress = F,
                   col_names = c("transcript", "gene", "gene_symbol","gene_full_symbol","chr","start","end","strand")) %>% 
    select(-c(transcript, gene_full_symbol, chr, start, end, strand)) %>% distinct()
  rownames(pbmc.kalisto.data) <- tr2g$gene_symbol[match(rownames(pbmc.kalisto.data), tr2g$gene)]
  kallisto_pbmc_before <- CreateSeuratObject(counts = pbmc.kalisto.data, project = "pbmc3k.kb", min.cells = 3, min.features = 100)
  ####
  starsolo_pbmc_before <- load_starsolo("/exports/sascstudent/thomas/output/STARsolo/pbmc3k/Solo.out/Gene/raw", "starsolo_pbmc")
  remove(pbmc.data, pbmc.kalisto.data, tr2g)
}

cellranger_pbmc_before
kallisto_pbmc_before
starsolo_pbmc_before

cellranger_pbmc_before[["percent.mt"]] <- PercentageFeatureSet(cellranger_pbmc_before, pattern = "^MT-")
kallisto_pbmc_before[["percent.mt"]] <- PercentageFeatureSet(kallisto_pbmc_before, pattern = "^MT-")
starsolo_pbmc_before[["percent.mt"]] <- PercentageFeatureSet(starsolo_pbmc_before, pattern = "^MT-")

Idents(cellranger_pbmc_before) <- "Cellranger"
Idents(kallisto_pbmc_before) <- "Kallisto"
Idents(starsolo_pbmc_before) <- "STARsolo"

p1 <- VlnPlot(cellranger_pbmc_before, features = "nFeature_RNA") + NoLegend() + ylim(0,3000) + xlab("") + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
p2 <- VlnPlot(kallisto_pbmc_before, features = "nFeature_RNA") + NoLegend() + ylim(0,3000) + xlab("") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
p3 <- VlnPlot(starsolo_pbmc_before, features = "nFeature_RNA") + NoLegend() + ylim(0,3000) + xlab("") + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
patch <- p1 | p2 | p3
patch

p1 <- VlnPlot(cellranger_pbmc_before, features = "percent.mt") + NoLegend() + ylim(0,60) + xlab("") + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
p2 <- VlnPlot(kallisto_pbmc_before, features = "percent.mt") + NoLegend() + ylim(0,60) + xlab("") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
p3 <- VlnPlot(starsolo_pbmc_before, features = "percent.mt") + NoLegend() + ylim(0,60) + xlab("") + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
patch <- p1 | p2 | p3
patch
remove(p1,p2,p3,patch)
remove(cellranger_pbmc_before,kallisto_pbmc_before,starsolo_pbmc_before)

#################################
cellranger_pbmc <- readRDS("/exports/sascstudent/thomas/R_output/cell_pbmc.rds")
kallisto_pbmc <- readRDS("/exports/sascstudent/thomas/R_output/kal_pbmc.rds")
starsolo_pbmc <- readRDS("/exports/sascstudent/thomas/R_output/starsolo_pbmc3k.rds")

cellranger_pbmc
kallisto_pbmc
starsolo_pbmc

VlnPlot(cellranger_pbmc, features = "nFeature_RNA", group.by = "orig.ident")
VlnPlot(kallisto_pbmc, features = "nFeature_RNA", group.by = "orig.ident")
VlnPlot(starsolo_pbmc, features = "nFeature_RNA", group.by = "orig.ident")

VariableFeaturePlot(cellranger_pbmc, pt.size = 0.2) #11714 non-variable
VariableFeaturePlot(kallisto_pbmc, pt.size = 0.2) #17008 non-variable
VariableFeaturePlot(starsolo_pbmc) #12902 non-variable

Project(cellranger_pbmc) <- "Cellranger"
Project(kallisto_pbmc) <- "Kallisto"
Project(starsolo_pbmc) <- "STARsolo"

cellranger_pbmc[["CR_seurat_clusters"]] <- cellranger_pbmc$seurat_clusters
kallisto_pbmc[["CR_seurat_clusters"]] <- recode_factor(kallisto_pbmc$seurat_clusters, "1"="0","2"="1","0"="2")
starsolo_pbmc[["CR_seurat_clusters"]] <- recode_factor(starsolo_pbmc$seurat_clusters, "2"="1","1"="2")

DimPlot(cellranger_pbmc, reduction = "umap", group.by = "CR_seurat_clusters", order = c(8,7,6,5,4,3,2,1,0))
DimPlot(kallisto_pbmc, reduction = "umap", group.by = "CR_seurat_clusters", order = c(8,7,6,5,4,3,2,1,0))
DimPlot(starsolo_pbmc, reduction = "umap", group.by = "CR_seurat_clusters", order = c(8,7,6,5,4,3,2,1,0))

#Rand index to compare the two clustering globally #############################
adj_rand_index_seurat <- function(project_1,project_2){
  project1.clusters <- as.data.frame(project_1[[]] %>% select(seurat_clusters))
  project2.clusters <- as.data.frame(project_2[[]] %>% select(seurat_clusters))
  rownames(project1.clusters) <- rownames(project1.clusters) %>% str_extract(. , pattern = ".+[ATCG]+")
  rownames(project2.clusters) <- rownames(project2.clusters) %>% str_extract(. , pattern = ".+[ATCG]+")
  com.cluster.barcodes <- sort(intersect(rownames(project1.clusters), rownames(project2.clusters)))
  project1_cluster <- as.integer(project1.clusters[com.cluster.barcodes,"seurat_clusters"])
  project2_cluster <- as.integer(project2.clusters[com.cluster.barcodes,"seurat_clusters"])
  
  print(paste0("The adjusted rand index of the clusters from ",Project(project_1)," and ",
               Project(project_2)," is ",adj.rand.index(project1_cluster, project2_cluster)))
}

rand_index_seurat <- function(project_1,project_2){
  project1.clusters <- as.data.frame(project_1[[]] %>% select(seurat_clusters))
  project2.clusters <- as.data.frame(project_2[[]] %>% select(seurat_clusters))
  rownames(project1.clusters) <- rownames(project1.clusters) %>% str_extract(. , pattern = ".+[ATCG]+")
  rownames(project2.clusters) <- rownames(project2.clusters) %>% str_extract(. , pattern = ".+[ATCG]+")
  com.cluster.barcodes <- sort(intersect(rownames(project1.clusters), rownames(project2.clusters)))
  
  project1_cluster <- as.integer(project1.clusters[com.cluster.barcodes,"seurat_clusters"])
  project2_cluster <- as.integer(project2.clusters[com.cluster.barcodes,"seurat_clusters"])
  
  print(paste0("The rand index of the clusters from ",Project(project_1)," and ",
               Project(project_2)," is ",rand.index(project1_cluster, project2_cluster)))
}

rand_index_seurat(cellranger_pbmc, kallisto_pbmc)
rand_index_seurat(cellranger_pbmc, starsolo_pbmc)
rand_index_seurat(kallisto_pbmc, starsolo_pbmc)

adj_rand_index_seurat(cellranger_pbmc, kallisto_pbmc)
adj_rand_index_seurat(cellranger_pbmc, starsolo_pbmc)
adj_rand_index_seurat(kallisto_pbmc, starsolo_pbmc)


#INHOUSE#####################################


kal_inhouse_before[["percent.mt"]] <- PercentageFeatureSet(kal_inhouse_before, pattern = "^MT-")
starsolo_inhouse_before[["percent.mt"]] <- PercentageFeatureSet(starsolo_inhouse_before, pattern = "^MT-")
cellranger_inhouse_before[["percent.mt"]] <- PercentageFeatureSet(cellranger_inhouse_before, pattern = "^MT-")

cellranger_inhouse_before[["method"]] <- "Cellranger"
kal_inhouse_before[["method"]] <- "Kallisto"
starsolo_inhouse_before[["method"]] <- "STARsolo"

cellranger_inhouse_before
kal_inhouse_before
starsolo_inhouse_before

Idents(cellranger_inhouse_before) <- "Cellranger"
Idents(kal_inhouse_before) <- "Kallisto"
Idents(starsolo_inhouse_before) <- "STARsolo"


p1 <- VlnPlot(cellranger_inhouse_before, features = "nFeature_RNA") + NoLegend() + ylim(0,10000) + xlab("") + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
p2 <- VlnPlot(kal_inhouse_before, features = "nFeature_RNA") + NoLegend() + ylim(0,10000) + xlab("") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
p3 <- VlnPlot(starsolo_inhouse_before, features = "nFeature_RNA") + NoLegend() + ylim(0,10000) + xlab("") + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
patch <- p1 | p2 | p3
patch

p1 <- VlnPlot(cellranger_inhouse_before, features = "percent.mt") + NoLegend() + ylim(0,100) + xlab("") + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
p2 <- VlnPlot(kal_inhouse_before, features = "percent.mt") + NoLegend() + ylim(0,100) + xlab("") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
p3 <- VlnPlot(starsolo_inhouse_before, features = "percent.mt") + NoLegend() + ylim(0,100) + xlab("") + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
patch <- p1 | p2 | p3
patch

{
  cellranger_inhouse <- readRDS("/exports/sascstudent/thomas/R_output/cellranger_inhouse.rds")
  kallisto_inhouse <- readRDS("/exports/sascstudent/thomas/R_output/hiPSC_kallisto_inhouse_mnn.rds")
  starsolo_inhouse <- readRDS("/exports/sascstudent/thomas/R_output/starsolo_inhouse.rds")
  original_paper <- readRDS("/exports/sascstudent/thomas/hiPSC_susana/hiPSC_complete_mnn.rds")
  
  Project(cellranger_inhouse) <- "cellranger"
  Project(kallisto_inhouse) <- "kallisto"
  Project(starsolo_inhouse) <- "starsolo"
  Project(original_paper) <- "publication"
  
  cellranger_inhouse[["CR_seurat_clusters"]] <- cellranger_inhouse$seurat_clusters
  kallisto_inhouse[["CR_seurat_clusters"]] <- kallisto_inhouse$seurat_clusters
  starsolo_inhouse[["CR_seurat_clusters"]] <- starsolo_inhouse$seurat_clusters
  original_paper[["CR_seurat_clusters"]] <- original_paper$seurat_clusters
}

rand_index_seurat(cellranger_inhouse, kallisto_inhouse)
rand_index_seurat(cellranger_inhouse, starsolo_inhouse)
rand_index_seurat(starsolo_inhouse, kallisto_inhouse)

adj_rand_index_seurat(cellranger_inhouse, kallisto_inhouse)
adj_rand_index_seurat(cellranger_inhouse, starsolo_inhouse)
adj_rand_index_seurat(starsolo_inhouse, kallisto_inhouse)



#Fix for cellranger suffix
get_cluster_count <- function(cell.cluster.num, project1, project2) {
  cell.cluster.barcodes <- project1[[]] %>%
    filter(CR_seurat_clusters == cell.cluster.num) %>%
    rownames() %>%
    str_extract(. , pattern = ".+[ATCG]+")
  project1.cluster.name <- paste0(Project(project1),".clust_",cell.cluster.num)
  test <- project2[[]]%>%
    select(CR_seurat_clusters) %>%
    rownames_to_column() %>%
    filter(rowname %in% cell.cluster.barcodes) %>%
    column_to_rownames() %>%
    group_by(CR_seurat_clusters) %>%
    summarise(project1.cluster.name = length(CR_seurat_clusters)) %>%
    column_to_rownames("CR_seurat_clusters")
  rownames(test) <- lapply(rownames(test), function(x) paste0(Project(project2), ".clust_", x))
  colnames(test) <- paste0(Project(project1),".clust_",cell.cluster.num)
  test <- test %>% rownames_to_column()
  return(test)
}
get_shared_cells <- function(project1, project2){
  clus_list <- lapply(levels(project1$CR_seurat_clusters), FUN = function(x) get_cluster_count(x, project1, project2))
  print(length(clus_list))
  print(levels(project1@meta.data$CR_seurat_clusters))
  clust_share <- data.frame(x=0:length(clus_list))
  rownames(clust_share) <- lapply(seq(1:(length(clus_list)+1))-1, function(x) paste0(Project(project2), ".clust_", x))
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

clust_share <- get_shared_cells(kallisto_inhouse, starsolo_inhouse)
clust_share <- clust_share[1:length(rownames(clust_share))-1,]

#Clusters similarity heatmap
{
  sample1.clusters <- colnames(clust_share)
  sample2.clusters <- rownames(clust_share)
  data <- expand_grid(sample1.clusters, sample2.clusters)
  for (i in seq(1:nrow(data))){
    data[i, "cells"] <- clust_share[as.character(data[i,"sample2.clusters"]), as.character(data[i,"sample1.clusters"])]
  }
  remove(sample1.clusters,sample2.clusters,i)
}

data %>%
  mutate(sample1.clusters=factor(sample1.clusters, levels = unique(sample1.clusters))) %>%
  mutate(sample2.clusters=factor(sample2.clusters, levels = unique(sample2.clusters))) %>%
  ggplot(aes(sample1.clusters, sample2.clusters, fill = cells)) + 
  geom_tile() + 
  geom_text(aes(label = cells), color = "#FFFFFF") +
  theme(axis.text.x = element_text(angle=70, hjust = 1))

#INHOUSE MARKERS###############################################################

cellranger_inhouse.markers <- FindAllMarkers(cellranger_inhouse, only.pos = T, 
                                             logfc.threshold = 0.25, min.pct = 0.25)
kallisto_inhouse.markers <- FindAllMarkers(kallisto_inhouse, only.pos = T, 
                                           logfc.threshold = 0.25, min.pct = 0.25)
starsolo_inhouse.markers <- FindAllMarkers(starsolo_inhouse, only.pos = T, 
                                           logfc.threshold = 0.25, min.pct = 0.25)
original_paper.markers <- FindAllMarkers(original_paper, only.pos = T, 
                                         logfc.threshold = 0.25, min.pct = 0.25)

gost_seurat <- function(markers, short_link = F){
  query_list <- list()
  clusters <- unique(markers$cluster)
  for (cl in clusters) {
    cluster <- paste0("cl_", cl)
    query_list[[cluster]] <- markers %>%
      filter(cluster == cl) %>%
      arrange(p_val_adj) %>%
      dplyr::select(gene)
  }
  summary(query_list)
  gostres <- gost(query_list, organism = "hsapiens",
       ordered_query = TRUE, user_threshold = 0.05, domain_scope = "annotated", 
       sources = c("GO:BP", "KEGG"), as_short_link = short_link)
  return(gostres)
}


##c('CO:MF', 'GO:CC', 'GO:BP', 'KEGG', 'REAC'))
cellranger_inhouse.gost <- gost_seurat(cellranger_inhouse.markers, F)
gost_seurat(cellranger_inhouse.markers, T)
kallisto_inhouse.gost <- gost_seurat(kallisto_inhouse.markers, F)
gost_seurat(kallisto_inhouse.markers, T)
starsolo_inhouse.gost <- gost_seurat(starsolo_inhouse.markers, F)
gost_seurat(starsolo_inhouse.markers, T)

#Filtered markers############################

cellranger_inhouse.markers.filt <- cellranger_inhouse.markers %>% 
  filter(p_val_adj < 0.01)
kallisto_inhouse.markers.filt <- kallisto_inhouse.markers %>% 
  filter(p_val_adj < 0.01)
starsolo_inhouse.markers.filt <- starsolo_inhouse.markers %>% 
  filter(p_val_adj < 0.01)
original_paper.markers.filt <- original_paper.markers %>%
  filter(p_val_adj < 0.01)

#Top markers#################################

cellranger_inhouse.markers.top <- cellranger_inhouse.markers.filt %>%
  group_by(cluster)
kallisto_inhouse.markers.top <- kallisto_inhouse.markers.filt %>%
  group_by(cluster)
starsolo_inhouse.markers.top <- starsolo_inhouse.markers.filt %>%
  group_by(cluster)
original_paper.markers.top <- original_paper.markers.filt %>%
  group_by(cluster)


# write.table(x = cellranger_inhouse.markers.top,
#             file = "./cellranger_topmarkers.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(x = kallisto_inhouse.markers.top,
#             file = "./kallisto_topmarkers.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(x = starsolo_inhouse.markers.top,
#             file = "./starsolo_topmarkers.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
#write.table(x = original_paper.markers.top,
#             file = "./original_paper_topmarkers.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


cellranger_inhouse.markers <- read.table(file = "./cellranger_topmarkers.tsv",
                                             sep = "\t", header = T)
kallisto_inhouse.markers <- read.table(file = "./kallisto_topmarkers.tsv",
                                           sep = "\t", header = T)
starsolo_inhouse.markers <- read.table(file = "./starsolo_topmarkers.tsv",
                                           sep = "\t", header = T)
original_paper.markers <- read.table(file = "./original_paper_topmarkers.tsv",
                                         sep = "\t", header = T)

{
  cellranger_inhouse.markers.top <- cellranger_inhouse.markers %>% 
    filter(p_val_adj<0.01) %>% group_by(cluster) %>% slice_min(n=50, order_by = p_val_adj)
  kallisto_inhouse.markers.top <- kallisto_inhouse.markers %>% 
    filter(p_val_adj<0.01) %>% group_by(cluster) %>% slice_min(n=50, order_by = p_val_adj)
  starsolo_inhouse.markers.top <- starsolo_inhouse.markers %>% 
    filter(p_val_adj<0.01) %>% group_by(cluster) %>% slice_min(n=50, order_by = p_val_adj)
  original_paper.markers.top <- original_paper.markers %>% 
    filter(p_val_adj<0.01) %>% group_by(cluster) %>% slice_min(n=50, order_by = p_val_adj)
}

##############################################################

marker_sim <- function(proj1markers, proj2markers) {
  total_weight <- 0
  static_weight <- 1/length(proj1markers$gene)
  for(i in seq(1:length(proj1markers$gene))) {
    if (proj1markers$gene[i] == proj2markers$gene[i]) {
      total_weight <- total_weight + static_weight
    } else {
      for (j in seq(1:length(proj2markers$gene))) {
        if (proj1markers$gene[i] == proj2markers$gene[j]) {
          dyn_weight <- static_weight * (1 - ((1 / length(proj2markers$gene) * abs(i - j))))
          total_weight <- total_weight + dyn_weight
        }
      }
    }
  }
  return(total_weight)
}

{
  #Custom similarity method
  sims <- vector(mode = "list", length = length(unique(kallisto_inhouse.markers.top$cluster)))
  for (nr in sort(unique(kallisto_inhouse.markers.top$cluster))) {
    proj1_markers <- kallisto_inhouse.markers.top %>% filter(p_val_adj<0.01) %>% 
      filter(cluster==nr) %>% slice_max(n = 50, order_by = p_val_adj)
    proj2_markers <- starsolo_inhouse.markers.top %>% filter(p_val_adj<0.01) %>% 
      filter(cluster==nr) %>% slice_head(n = 50)
    sims[(nr+1)] <- marker_sim(proj1_markers, proj2_markers)
  }
  sims
}

jac_index <- function(a, b){
  res <- vector(mode = "list", length = length(unique(a$cluster)))
  for (nr in sort(unique(a$cluster))) {
    a_cl <- filter(a, cluster==nr)
    b_cl <- filter(b, cluster==nr)
    inter_size <- length(intersect(a_cl$gene, b_cl$gene))
    union_size <- length(a_cl$gene) + length(b_cl$gene) - inter_size
    res[[nr+1]] <- inter_size / union_size
  }
  return(res)
}

{
  #Jaccard index
  paste(jac_index(cellranger_inhouse.markers.top, kallisto_inhouse.markers.top), collapse = ", ")
  paste(jac_index(cellranger_inhouse.markers.top, starsolo_inhouse.markers.top), collapse = ", ")
  paste(jac_index(cellranger_inhouse.markers.top, original_paper.markers.top), collapse = ", ")
  
  paste(jac_index(kallisto_inhouse.markers.top, starsolo_inhouse.markers.top), collapse = ", ")
}


{
  shared_barcodes <- intersect(Cells(cellranger_inhouse), Cells(original_paper))
  barcode_counts <- merge(cellranger_inhouse[[c("nCount_RNA","nFeature_RNA","sample_name","orig.ident")]], 
                          y = original_paper[[c("nCount_RNA","nFeature_RNA","sample_name","run")]], 
                          by="row.names", all=T)
  colnames(barcode_counts) <- c("cells", "nCount_RNA.rerun","nFeature_RNA.rerun", "samples.rerun","timepoint.rerun",
                                "nCount_RNA.orig","nFeature_RNA.orig","samples.orig","timepoint.orig")
  
  barcode_counts["nCount.diff"] <- barcode_counts["nCount_RNA.rerun"] - barcode_counts["nCount_RNA.orig"]
  barcode_counts["nFeature.diff"] <- barcode_counts["nFeature_RNA.rerun"] - barcode_counts["nFeature_RNA.orig"]
  ggplot(barcode_counts, aes(y=nCount.diff, fill=timepoint.rerun)) + 
    geom_boxplot()
  ggplot(barcode_counts, aes(y=nFeature.diff, fill=timepoint.rerun)) + 
    geom_boxplot()
  
  lowercount_cells <- select(barcode_counts, c(cells, nCount.diff)) %>% filter(nCount.diff < -8) %>% select(cells)
  DimPlot(cellranger_inhouse, cells.highlight = lowercount_cells$cells)
  DimPlot(original_paper, cells.highlight = lowercount_cells$cells)
  summary(barcode_counts$nFeature.diff)
  
  lowerfeature_cells <- select(barcode_counts, c(cells, nFeature.diff)) %>% filter(nFeature.diff < -8) %>% select(cells)
  DimPlot(cellranger_inhouse, cells.highlight = lowerfeature_cells$cells)
  DimPlot(original_paper, cells.highlight = lowerfeature_cells$cells)
}






