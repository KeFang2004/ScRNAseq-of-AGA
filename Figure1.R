HFC <- obj[, Idents(obj) %in% "HFC"]

HFC <- NormalizeData(HFC, verbose = F) 
HFC <- FindVariableFeatures(HFC, selection.method = 'vst', nfeatures = 2000)
HFC <- ScaleData(HFC, vars.to.regress = "percent.mt")
HFC <- RunPCA(HFC, features = VariableFeatures(object = HFC)) 
HFC <- RunHarmony(HFC,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
HFC <- RunUMAP(HFC, reduction = "harmony", dims = 1:30,reduction.name = "umap")
HFC <- FindNeighbors(HFC, dims = 1:30)
#HFC_copy <- HFC
#HFC_copy <- FindClusters(HFC_copy, resolution = seq(0.1,1.3,0.2))
#library(clustree)
#clustree(HFC_copy@meta.data, prefix = "RNA_snn_res.")
HFC <- FindClusters(HFC, resolution = 0.4)
# Look at cluster IDs of the first 5 cells
head(Idents(HFC), 5)
table(HFC$seurat_clusters) 
DimPlot(HFC, reduction = 'umap', label = TRUE)

new.cluster.ids <- c("Outer layer cluster", "Inner layer cluster", "Outer layer cluster", "Inner layer cluster", "Outer layer cluster", "Outer layer cluster",
                     "Outer layer cluster", "Inner layer cluster", "Outer layer cluster", "Inner layer cluster", "Outer layer cluster","Outer layer cluster",
                     "Outer layer cluster","Inner layer cluster","Outer layer cluster","Outer layer cluster",
                     "Inner layer cluster","Outer layer cluster")
names(new.cluster.ids) <- levels(HFC)
HFC <- RenameIdents(HFC, new.cluster.ids)
current.cluster.ids <- c(0:17)
#替换分群序号为注释
HFC@meta.data$cell_type = plyr::mapvalues(x = HFC@meta.data[,"seurat_clusters"], from = current.cluster.ids, to = new.cluster.ids)
clusterCols <- c("#ed4437","#1f78b4")
names(clusterCols) <- c("Inner layer cluster", "Outer layer cluster")
DimPlot(HFC, reduction = "umap", cols = clusterCols, label = TRUE, pt.size = 0.3)+
  coord_fixed(ratio = 1)
HFC.markers <- FindAllMarkers(HFC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 2.5) 
library(dplyr)
HFC.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) 
top10 <- HFC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
FeaturePlot(HFC, features = "KRT25", reduction = "umap", cols = c("#F8F8F8","#ed4437"),pt.size = 0.3,order = T)+
  coord_fixed(ratio = 1)


library(UCell)
markers <- list()
markers$sbORS <- c("KRT5","KRT6A","KRT6B","KRT6C","KRT17","MSX2-","KRT15-","ITGA6-","KRT35-","KRT85-")
markers$bORS <- c("KRT5","ITGA6","KRT15","MSX2-","KRT6A-","KRT6B-","KRT6C-","KRT17-","KRT35-","KRT85-")
markers$CP <- c("KRT75","CST6","KRT79","KRT6A-","KRT6B-","KRT6C-","KRT17-","MSX2-","KRT15-","ITGA6-","KRT35-","KRT85-")
markers$IRS <- c("MSX2","KRT73","KRT71","KRT25","KRT28","KRT5-","KRT85-")
markers$Cuticle <- c("KRT82")
markers$Cortex <- c("MSX2","KRT85","KRT5-","KRT73-","KRT71-","KRT82-","KRT14-")
marker_score <- AddModuleScore_UCell(HFC, features=markers)
library(stringr)
library(ggplot2)
library(viridis)
a <- colnames(marker_score@meta.data) %>% str_subset("_UCell")
options(repr.plot.width=8, repr.plot.height=6)
FeaturePlot(marker_score, reduction = "umap", features = a,
            ncol = 2, order = T,
            cols = c("#F8F8F8","#077431"), pt.size = 0.3)
write.csv(HFC.markers, "markers_OL_IL.csv")
