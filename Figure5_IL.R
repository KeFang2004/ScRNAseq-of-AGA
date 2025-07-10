IL <- HFC[, Idents(HFC) %in% "Inner layer cluster"]
IL <- NormalizeData(IL, verbose = F) 
IL <- FindVariableFeatures(IL, selection.method = 'vst', nfeatures = 2000)
IL <- ScaleData(IL, vars.to.regress = "percent.mt")
IL <- RunPCA(IL, features = VariableFeatures(object = IL)) 
IL <- RunHarmony(IL,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
IL <- RunUMAP(IL, reduction = "harmony", dims = 1:30,reduction.name = "umap")
IL <- FindNeighbors(IL, dims = 1:30)
#IL_copy <- IL
#IL_copy <- FindClusters(IL_copy, resolution = seq(0.1,1.3,0.2))
#library(clustree)
#clustree(IL_copy@meta.data, prefix = "RNA_snn_res.")
IL <- FindClusters(IL, resolution = 0.3)
DimPlot(IL, reduction = 'umap', label = TRUE)
genes_to_check = c("MSX2","KRT85", "KRT35", "KRT36", "MKI67", "LEF1", "KRT82", "KRT73","GATA3","TCHH")
DotPlot(IL,group.by = 'seurat_clusters', features = genes_to_check,cluster.idents = T) +
  scale_color_gradientn(colours = c('#64cccf', 'white','#ea5c6f'))+
  theme(axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))
DotPlot(IL,group.by = 'seurat_clusters', features = genes_to_check,cluster.idents = T) +
  scale_color_gradientn(colours = c('#FFCC33', '#66CC66','#336699','#330066'))+
  theme(axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))
il.markers <- FindAllMarkers(IL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1.5) 
library(dplyr)
il.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) 
top10 <- il.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
FeaturePlot(IL, features = "SOX21", reduction = "umap", pt.size = 0.5)

new.cluster.ids <- c("CO1", "HMC", "CO2", "IRS", "MIC", "CU",
                     "CO4", "CO3", "MED")
names(new.cluster.ids) <- levels(IL)
IL <- RenameIdents(IL, new.cluster.ids)
current.cluster.ids <- c(0,1,2,3,4,5,6,7,9)
#替换分群序号为注释
IL@meta.data$cell_type = plyr::mapvalues(x = IL@meta.data[,"seurat_clusters"], from = current.cluster.ids, to = new.cluster.ids)
clusterCols <- c("#d4a6a8","#af93c4","#8660a8","#6a3d9a","#b81316","#5c9e43","#f48521","#f9b769","#c6b598")
names(clusterCols) <- c("CO1",
                        "CO2",
                        "CO3",
                        "CO4",
                        "MED",
                        "CU",
                        "IRS",
                        "MIC",
                        "HMC")
DimPlot(IL, reduction = "umap", label = TRUE, cols = clusterCols, pt.size = 0.3)+
  coord_fixed(ratio = 1)
Idents(IL) <- factor(Idents(IL), levels =  c("CO1",
                                             "CO2",
                                             "CO3",
                                             "CO4",
                                             "MED",
                                             "CU",
                                             "IRS",
                                             "MIC",
                                             "HMC") )
IL1 <- IL
IL1 <- RenameIdents(IL1, CO1 = "CO", CO2 = "CO", CO3 = "CO", CO4 = "CO")
il1.markers <- FindAllMarkers(IL1, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 2.0) 
top8 <- il1.markers %>% group_by(cluster) %>% top_n(n = 8, wt = avg_log2FC)

library(Seurat)
library(ggplot2)
modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  p <- VlnPlot(obj, features = feature, pt.size = pt.size, ...) +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin) +
    scale_y_continuous(labels = scales::comma)  # 添加纵轴标签格式化
  return(p)
}

## main function
StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj, feature = x, ...))
  plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
    theme(axis.text.x = element_text(), axis.ticks.x = element_line())
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

StackedVlnPlot(IL, c("XXbac-BPG32J3.19",'EFHD1', 'KRTAP11-1', 'SELENBP1', 'DAPL1', 'KRT31','KRT85',"DSG4"), pt.size=0, cols=my36colors)
#配色方案
my36colors <- c("#d4a6a8","#af93c4","#8660a8","#6a3d9a","#b81316","#5c9e43","#f48521","#f9b769","#c6b598")


library(destiny)
set.seed(100)
# 进行扩散映射
dm = destiny::DiffusionMap(data = t(as.matrix(IL@assays$RNA$scale.data[VariableFeatures(IL),])) ,n_pcs = 30)
# 绘制扩散映射的前两个主成分
plot(dm, 1:2)
# 指定HMC作为起始细胞
tips = rownames(IL@meta.data)[IL@meta.data$cell_type == "HMC"][1]
# 获取起始细胞在扩散映射中的索引
tips_index <- which(rownames(dm@eigenvectors) == tips)
plot(eigenvalues(dm), ylim = 0:1, pch = 20, xlab = 'Diffusion component (DC)', ylab = 'Eigenvalue')
# 计算伪时间
dpt <- DPT(dm, tips = tips_index)
# 绘制伪时间轨迹
plot(dpt)
tmp <- data.frame(DC1 = dm$DC1,
                  DC2 = dm$DC2,
                  Timepoint = IL$cell_type,
                  dpt = dpt$DPT1)
library(ggthemes)
library(ggbeeswarm)
# 定义颜色
my_colors <- c("CO1" = "#d4a6a8", "CO2" = "#af93c4", "CO3" = "#8660a8", "CO4" = "#6a3d9a",
               "MED" = "#b81316", "CU" = "#5c9e43", "IRS" = "#f48521", "MIC" = "#f9b769", "HMC" = "#c6b598")
ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint)) +  
  geom_point(alpha = 0.7) +
  scale_color_tableau() + 
  scale_color_manual(values = my_colors) +
  xlab("Diffusion component 1") +   
  ylab("Diffusion component 2") +  
  theme_classic()

# Pseudotime analysis
library(monocle3)
mat <- GetAssayData(IL, slot = "data")  
cellInfo <- IL@meta.data  
geneInfo <- data.frame(gene_short_name = rownames(mat), row.names = rownames(mat))  
cds <- new_cell_data_set(expression_data = mat,  
                         cell_metadata = cellInfo,  
                         gene_metadata = geneInfo)  
cds <- preprocess_cds(cds, num_dim = 50)
reducedDim(cds, "UMAP") <- Embeddings(IL, "umap")
cds = cluster_cells(cds, cluster_method = 'louvain')
cds = learn_graph(cds, use_partition=T, verbose=T, learn_graph_control=list(
  minimal_branch_len=10
)) #Trajectory learning
# 获取HMC细胞的索引
hmc_cell_ids <- rownames(colData(cds))[colData(cds)$cell_type == "HMC"]
# 找到HMC细胞在UMAP图中的索引
closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[hmc_cell_ids, ]))))]
# 将起点指定为HMC细胞
cds <- order_cells(cds, root_pr_nodes = root_pr_nodes)
# 使用Seurat的UMAP信息，这样可以与Seurat对象的细胞分布保持一致
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(IL, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
##使用Seurat的UMAP信息，这样可以与Seurat对象的细胞分布保持一致
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(IL, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
clusterCols <- c("#d4a6a8","#af93c4","#8660a8","#6a3d9a","#b81316","#5c9e43","#f48521","#f9b769","#c6b598")
names(clusterCols) <- c("CO1", "CO2", "CO3", 
                        "CO4", "MED", "CU", 
                        "IRS", "MIC", "HMC")
plot_cells(cds, color_cells_by = "pseudotime",
           show_trajectory_graph=T,cell_size = 0.55) + plot_cells(cds,
                                                 color_cells_by = "cell_type",
                                                 label_cell_groups=FALSE,
                                                 label_leaves=FALSE,
                                                 label_branch_points=FALSE,
                                                 graph_label_size=1,
                                                 cell_size = 0.55)+
  scale_color_manual(values = clusterCols) +
  scale_fill_manual(values = clusterCols)

library(loonR)
modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores =6)
genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.25))
mat <- pre_pseudotime_matrix(cds, genes)
# check
head(mat[1:5,1:5])

library(ClusterGVis)
mat <- as.data.frame(mat)
# kmeans
ck <- clusterData(mat,
                  cluster.method = "kmeans",
                  cluster.num = 9)

# add line annotation
library(org.Hs.eg.db)
enrich <- enrichCluster(object = ck,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        pvalueCutoff = 0.05,
                        topn = 5,
                        seed = 5201314)
pdf('Figure5_monocle3_heatmap.pdf',height = 10,width = 8,onefile = F)
visCluster(object = ck,
           plot.type = "both",
           add.sampleanno = F,
           markGenes = sample(rownames(mat),50,replace = F),
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left")
visCluster(object = ck,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = sample(rownames(mat),50,replace = F),
           markGenes.side = "left",
           annoTerm.data = enrich,
           add.bar = T,
           line.side = "left")
dev.off()

install.packages("jjAnno")
