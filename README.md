# Analyzing transcriptomic changes upon stem cell to neuron differentiation in the developing cerebral cortex
Pipeline to analyze single-cell data from Seurat and perform trajectory analysis with Monocle3
```
library(Seurat)
library(tidyverse)
library(dplyr)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(ggridges)
```
## Read 10X dataset(feature barcode matrix) and create a Seurat Object
```
data <- Read10X(data.dir = "D:/ChIP seq and RNA seq analysis/Arlotta_lab scRNA/E17.5.filtered_feature_bc_matrix")
data <- CreateSeuratObject(counts = data, project='E17.5', min.cells = 3, min.features = 200)
```
## Check and Visualize QC metrics like mitochonrial RNA reads.
```
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```
![image](https://github.com/mahibose/Seurat_to_Monocle3_v2/assets/60874215/b7af7b9f-ddc2-45cd-b094-229407caf1f0)

## Visualizing the quality of reads per cell
```
FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth()
```
![image](https://github.com/mahibose/Seurat_to_Monocle3_v2/assets/60874215/55550268-bed6-4d23-bf19-ef87a408395f)

## Filtering and removing low-quality cells
```
data <- subset(data, subset = nFeature_RNA >200 & nFeature_RNA <2500 & percent.mt <5)
```
## Normalizing data
```
data <- NormalizeData(data)
```
## Identify highly variable features
```
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(data), 10)
top10
```
## Plot variable features
```
plot1 <- VariableFeaturePlot(data)
LabelPoints(plot = plot1, points = top10, repel = T)
```
![image](https://github.com/mahibose/Seurat_to_Monocle3_v2/assets/60874215/fecdcaf3-38b9-41e1-bc2a-2171f308c1df)

## Scaling data
```
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
```
## Perform dimensionality reduction
```
data <- RunPCA(data, features = VariableFeatures(object = data
DimHeatmap(data, dims = 1:6, cells = 500, balanced = T)
```
![image](https://github.com/mahibose/Seurat_to_Monocle3_v2/assets/60874215/0e6e160c-87b9-48fb-9a84-76fe8d1a4d33)

## Determine the Dimensionality of the data
```
ElbowPlot(data)
```
![image](https://github.com/mahibose/Seurat_to_Monocle3_v2/assets/60874215/6fc31c27-d244-46d0-85f9-396aabd73014)

## Clustering data
```
data <- FindNeighbors(data, dims = 1:20)
```
## Understanding the resolution of the clustering
```
data <- FindClusters(data, resolution = c(0.3, 0.5, 0.7, 1))
head(data@meta.data)
p1 <- DimPlot(data, group.by = "RNA_snn_res.0.3", label = T)
p2 <- DimPlot(data, group.by = "RNA_snn_res.0.5", label = T)
p3 <- DimPlot(data, group.by = "RNA_snn_res.0.7", label = T)
p4 <- DimPlot(data, group.by = "RNA_snn_res.1", label = T)
p1 + p2 + p3 + p4
```
![image](https://github.com/mahibose/Seurat_to_Monocle3_v2/assets/60874215/623f8849-f089-44e3-8f8c-b9429aa03627)

## Setting Identity of Clusters
```
Idents(data) <- "RNA_snn_res.0.7"
```
## Non-linear Dimensionality Reduction: UMAP
```
data <- RunUMAP(data, dims = 1:20)
DimPlot(data, reduction = "umap")
```
![image](https://github.com/mahibose/Seurat_to_Monocle3_v2/assets/60874215/f5e613ba-9309-42df-9542-ca6ae890c687)

## Looking at Variable Features to identify cluster markers
```
VlnPlot(data, features = c("Pax6", "Rbfox1"), slot = "counts", log = TRUE)
```
![image](https://github.com/mahibose/Seurat_to_Monocle3_v2/assets/60874215/9fb1b314-f5a2-490b-b236-89c6bacb7e9f)

## Looking at cluster markers to identify UMAP clusters
```
FeaturePlot(data, features = c("Pax6",  "Eomes", "Aldh1l1",
                               "Tbr1",  "Olig2", "Sox2", "Cux2", "Neurog2"))
```
![image](https://github.com/mahibose/Seurat_to_Monocle3_v2/assets/60874215/d5e1a256-a2d1-4592-b7ae-b9b8f869e15a)

## Rename Clusters
```
new.cluster.ids <- c("A progenitors","B progenitors", "Cnn2+ Msx1+ population" , "B progenitors", "Putative early neuronal population", "Lmx1a+ population",
                      "Intermediate Projenitors")
names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5)
```
# Trajectory Analysis with Monocle3
## Converting Seurat object to cell dataset object for Monocle3

```
cds <- as.cell_data_set(data)
```
## Get cell metadata
```
head(colData(cds))
```
## Get feature/gene metadata
```
fData(CDs)
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name <- rownames(fData(cds))
head(fData(cds))
```
## Get counts
```
head(counts(CDs))
```
## Retrieve clustering information from the Seurat object
### 1. Assign partitions
```
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
```
### 2. Assign cluster information
```
list.cluster <- data@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
```
### 3. Assign UMAP coordinates
```
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- data@reductions$umap@cell.embeddings
```
## Plot
```
cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, 
           group_label_size = 5) + theme(legend.position = "right")
cluster.before.traj
```
![image](https://github.com/mahibose/Seurat_to_Monocle3_v2/assets/60874215/4cde983e-bcd0-4025-9b1c-cc4dae2c871c)

## Learn Trajectory
```
cds <- learn_graph(cds, use_partition = F)
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,
           label_branch_points = T, label_roots = T, label_leaves = F,
           group_label_size = 5)

```
![image](https://github.com/mahibose/Seurat_to_Monocle3_v2/assets/60874215/1cce7924-44c6-4b7a-bdcf-168f1c0d7d31)

## Order cells in Pseudotime
```
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, clusters(cds) == 5]))
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)
```
![image](https://github.com/mahibose/Seurat_to_Monocle3_v2/assets/60874215/a9eb6e67-934b-4caa-93e5-554096f1d9eb)

## Cells ordered by Monocle3 Pseudotime
```
head(pseudotime(cds), 10)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, seurat_clusters, fill = seurat_clusters)) + geom_boxplot()
```
![image](https://github.com/mahibose/Seurat_to_Monocle3_v2/assets/60874215/02b4cbda-5968-42c2-9003-fdcaa93bcaa7)

```
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) + geom_boxplot()
```
![image](https://github.com/mahibose/Seurat_to_Monocle3_v2/assets/60874215/49d8427b-a483-4ab8-bef0-f709c8188742)

## Find genes that change as a function of pseudotime
```
deg <- graph_test(cds, neighbor_graph = "principal_graph")
deg %>% arrange(q_value) %>% filter(status == "OK") %>% head()
FeaturePlot(data, features = c("Rgs20", "Bend6", "Mcm3", "B3gat2"))
```
![image](https://github.com/mahibose/Seurat_to_Monocle3_v2/assets/60874215/015dc589-5aed-42f6-b0c8-94b28c9b1612)

## Add pseudotime values into the SeuratObject
```
data$pseudotime <- pseudotime(cds)
FeaturePlot(data, features = "pseudotime")
```
![image](https://github.com/mahibose/Seurat_to_Monocle3_v2/assets/60874215/f0df7dfc-ff43-4300-aa75-d1b8be437b73)
```
RidgePlot(data, features = c("Sox2", "Eomes", "Neurod2"), sort = T, idents = c("5", "6", "0", "1", "7"))
```
![image](https://github.com/mahibose/Seurat_to_Monocle3_v2/assets/60874215/80c1b5b8-66f1-4336-8b1e-b0d0e04a71d8)
```
my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("Sox2", "Eomes", "Neurod2"))) 
cds_subset <- cds[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime" )
```
![image](https://github.com/mahibose/Seurat_to_Monocle3_v2/assets/60874215/3aa5fb62-f6c6-4b1f-a788-3e17895424cb)
