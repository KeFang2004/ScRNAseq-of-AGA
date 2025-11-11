HMC <- IL[, Idents(IL) %in% 'HMC']
HMC <- NormalizeData(HMC, verbose = F) 
HMC <- FindVariableFeatures(HMC, selection.method = 'vst', nfeatures = 2000)
HMC <- ScaleData(HMC, vars.to.regress = "percent.mt")
HMC <- RunPCA(HMC, features = VariableFeatures(object = HMC)) 
HMC <- RunHarmony(HMC,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
HMC <- RunUMAP(HMC, reduction = "harmony", dims = 1:30,reduction.name = "umap")
HMC <- FindNeighbors(HMC, dims = 1:30)
HMC_copy <- HMC
HMC_copy <- FindClusters(HMC_copy, resolution = seq(0.1,0.9,0.2))
library(clustree)
clustree(HMC_copy@meta.data, prefix = "RNA_snn_res.")
HMC <- FindClusters(HMC, resolution = 0.4)
HMC.markers <- FindAllMarkers(HMC, only.pos = TRUE, min.pct = 0.4, logfc.threshold = 1.0) 
DimPlot(HMC, reduction = 'umap', label = TRUE)

clusterCols <- c("#ea5c6f", "#f7905a", "#e187cb","#fb948d")
names(clusterCols) <- c(0,1,2,3)
DimPlot(HMC, reduction = "umap", label = F, cols = clusterCols, pt.size = 1)+
  coord_fixed(ratio = 1)



FeaturePlot(HMC, features = "SERPINH1", reduction = "umap",cols = c('#daddf5','#ea5c6f'),pt.size = 1,order = T)+
  coord_fixed(ratio = 1)
library(UCell)
markers <- list()
markers$CO <- c("EFHD1","SELENBP1","DAPL1","KRT31","KRT85","DSG4","KRT82-","CRAT-","SLC27A6-","S100A3-","KRT32-"
                ,"KRT28-","KRT71-","CTSC-","KRT27-","TCHH-")
markers$MED <- c("ALPL","GPRC5A","TSPAN1","CTSH","TSPAN15","ALDH1A3","ZFHX3","KRT75","DSG4-")
markers$CU <- c("KRT82","CRAT","SLC27A6","KRT32","S100A3","TGM1","GNG11","PROCR","EFHD1-")
markers$IRS <- c("KRT28","KRT71","CTSC","KRT27","TCHH","KRT73","DEGS2","KRT25","DSG4-")
marker_score <- AddModuleScore_UCell(HMC, features=markers)
library(stringr)
library(ggplot2)
library(viridis)
a <- colnames(marker_score@meta.data) %>% str_subset("_UCell")
options(repr.plot.width=8, repr.plot.height=6)
FeaturePlot(marker_score, reduction = "umap", features = a,
            ncol = 2, order = T,
            cols = c('#F7F7F7','#a13037'), pt.size = 0.3)
write.csv(HFC.markers, "markers_OL_IL.csv")






library(destiny)
set.seed(100)
# 进行扩散映射
dm = destiny::DiffusionMap(data = t(as.matrix(HMC@assays$RNA$scale.data[VariableFeatures(HMC),])) ,n_pcs = 30)
# 绘制扩散映射的前两个主成分
plot(dm, 1:2)
# 指定起始细胞
tips = rownames(HMC@meta.data)[HMC@meta.data$seurat_clusters == 3][1]
# 获取起始细胞在扩散映射中的索引
tips_index <- which(rownames(dm@eigenvectors) == tips)
plot(eigenvalues(dm), ylim = 0:1, pch = 20, xlab = 'Diffusion component (DC)', ylab = 'Eigenvalue')
# 计算伪时间
dpt <- DPT(dm, tips = tips_index)
# 绘制伪时间轨迹
plot(dpt)
tmp <- data.frame(DC1 = dm$DC1,
                  DC2 = dm$DC2,
                  DC3 = dm$DC3,
                  Timepoint = HMC$seurat_clusters,
                  dpt = dpt$DPT1)
library(ggthemes)
library(ggbeeswarm)
my_colors <- c("#ea5c6f", "#f7905a", "#e187cb","#fb948d")
ggplot(tmp, aes(x = DC1, y = DC2, z = DC3, colour = Timepoint)) +  
  geom_point() + scale_color_tableau() +   
  xlab("Diffusion component 1") +   ylab("Diffusion component 2") + zlab("Diffusion component 3") + theme_classic()
library(plotly)
plot_ly(tmp, x = ~DC1, y = ~DC2, z = ~DC3, color = ~Timepoint, colors = my_colors) %>%
  add_markers() %>%
  layout(scene = list(
    xaxis = list(title = 'Diffusion component 1'),
    yaxis = list(title = 'Diffusion component 2'),
    zaxis = list(title = 'Diffusion component 3')
  ))
pdf("diffusion_map_3d.pdf", width = 8, height = 6)

# 绘制三维散点图
with(tmp, plot3d(DC1, DC2, DC3, col = my_colors[Timepoint], size = 2, pch = 19, xlab = "DC1", ylab = "DC2", zlab = "DC3"))

# 设置视图
rgl::view3d(theta = 30, phi = 30)  # 调整视角



# Pseudotime analysis
library(monocle3)
mat <- GetAssayData(HMC, slot = "data")  
cellInfo <- HMC@meta.data  
geneInfo <- data.frame(gene_short_name = rownames(mat), row.names = rownames(mat))  
cds <- new_cell_data_set(expression_data = mat,  
                         cell_metadata = cellInfo,  
                         gene_metadata = geneInfo)  
cds <- preprocess_cds(cds, num_dim = 50)
reducedDim(cds, "UMAP") <- Embeddings(HMC, "umap")
cds = cluster_cells(cds, cluster_method = 'louvain')
cds = learn_graph(cds, use_partition=T, verbose=T, learn_graph_control=list(
  minimal_branch_len = 4
)) #Trajectory learning

start = c(3)
closest_vertex = cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex = as.matrix(closest_vertex[colnames(cds), ])
root_pr_nodes = igraph::V(principal_graph(cds)[["UMAP"]])$name
cell_ids <- rownames(colData(cds),"seurat_clusters" == 3)
flag = closest_vertex[cell_ids,]
flag = as.numeric(names(which.max(table(flag))))
root_pr_nodes = root_pr_nodes[flag]
cds = order_cells(cds, root_pr_nodes=root_pr_nodes)
##使用Seurat的UMAP信息，这样可以与Seurat对象的细胞分布保持一致
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(HMC, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
clusterCols <- c("#ea5c6f", "#f7905a", "#e187cb","#fb948d")
names(clusterCols) <- c(0,1,2,3)
plot_cells(cds, color_cells_by = "pseudotime",
           show_trajectory_graph=T,
           cell_size = 1,) + plot_cells(cds,
                                                 color_cells_by = "seurat_clusters",
                                                 label_cell_groups=FALSE,
                                                 label_leaves=FALSE,
                                                 label_branch_points=FALSE,
                                                 cell_size = 1,
                                                 graph_label_size=1)+
  scale_color_manual(values = clusterCols) +
  scale_fill_manual(values = clusterCols)

cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
get_earliest_principal_node <- function(cds, time_bin= "3"){
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == time_bin)
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
par(mfrow=c(1,1))
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds_3d))
clusterCols <- c("#ea5c6f", "#f7905a", "#e187cb","#fb948d")
names(clusterCols) <- c(0,1,2,3)
plot_cells_3d(cds_3d, color_cells_by="seurat_clusters",color_palette = clusterCols,cell_size = 50)
cds_3d_plot_obj
pdf(file = "Figure5_3D_plot.pdf", width = 6, height = 4)
