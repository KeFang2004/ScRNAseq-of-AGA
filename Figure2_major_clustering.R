FeaturePlot(obj, features = "COL17A1", reduction = "umap", pt.size = 0.1)
new.cluster.ids <- c("EPI", "EPI", "HFC", "EPI", "HFC", "EPI", "HFC", "EPI", "EPI", "EPI", "HFC","HFC",
                     "HFC","ML","HFC","HFC","IMM","VASC","FIB","IMM","HFC","EPI","HFC","MISC","EPI")
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
current.cluster.ids <- c(0:24)
#替换分群序号为注释
obj@meta.data$cell_type = plyr::mapvalues(x = obj@meta.data[,"seurat_clusters"], from = current.cluster.ids, to = new.cluster.ids)
clusterCols <- c("#6B5382", "#C1719C", "#D0ECD7",
                 "#A3E0BB", "#367FA4", "#D9BAD0", 
                 "#4FC3B0")
names(clusterCols) <- c("ML", "IMM", "EPI", "HFC", "MISC", "VASC", "FIB")
DimPlot(obj, reduction = "umap", cols = clusterCols, label = TRUE, pt.size = 0.1)

clusterCols <- c("#ec8447", "#9cd9c0")
names(clusterCols) <- c("Forehead", "Occipital")
DimPlot(obj, reduction = "umap",label = TRUE, cols = clusterCols, split.by = "group", pt.size = 0.01)
Idents(obj) <- obj@meta.data$group
