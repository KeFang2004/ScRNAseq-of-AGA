OL <- HFC[, Idents(HFC) %in% "Outer layer cluster"]
OL=CreateSeuratObject(counts = OL@assays$RNA@layers$counts,
                       meta.data = OL@meta.data)
OL <- NormalizeData(OL, verbose = F) 
OL <- FindVariableFeatures(OL, selection.method = 'vst', nfeatures = 2000)
OL <- ScaleData(OL, vars.to.regress = "percent.mt")
OL <- RunPCA(OL, features = VariableFeatures(object = OL)) 
OL <- RunHarmony(OL,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
OL <- RunUMAP(OL, reduction = "harmony", dims = 1:30,reduction.name = "umap")
OL <- FindNeighbors(OL, dims = 1:30)
OL <- OL %>% RunTSNE(reduction = "harmony", dims = 1:30)
DimPlot(OL, reduction = "umap", label = TRUE, pt.size = 0.3)+
  coord_fixed(ratio = 1)
#OL_copy <- OL
#OL_copy <- FindClusters(OL_copy, resolution = seq(0.1,0.5,0.1))
#library(clustree)
#clustree(OL_copy@meta.data, prefix = "RNA_snn_res.")
OL <- FindClusters(OL, resolution = 0.3)

FeaturePlot(OL, features = "RUFY3", reduction = "umap",pt.size = 0.3,order = T)

ol.markers <- FindAllMarkers(OL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1.5) 
library(dplyr)
ol.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) 
top10 <- ol.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


new.cluster.ids <- c("sbORS1", "BG", "bORS3", "sbORS2", "bORS1", "bTAC",
                     "bORS3", "sbORS1", "CP", "bORS2", "bORS1","sbORS1","bORS4")
names(new.cluster.ids) <- levels(OL)
OL <- RenameIdents(OL, new.cluster.ids)
current.cluster.ids <- c(0:12)
#替换分群序号为注释
OL@meta.data$cell_type = plyr::mapvalues(x = OL@meta.data[,"seurat_clusters"], from = current.cluster.ids, to = new.cluster.ids)
DimPlot(OL, reduction = "umap", label = TRUE, pt.size = 0.3)+
  coord_fixed(ratio = 1)

ORS <- OL[, Idents(OL) %in% c("sbORS1","sbORS2","bORS1","bORS2","bORS3","bORS4")]
DimPlot(ORS, reduction = "umap", label = TRUE, pt.size = 0.3,split.by = "group")+
  coord_fixed(ratio = 1)
Idents(ORS) <- factor(Idents(ORS), levels =  c("sbORS1","sbORS2","bORS1","bORS2","bORS3","bORS4"))


###   GO analysis   ###
library(clusterProfiler)
library(org.Hs.eg.db)
ors.markers <- ors.markers[ors.markers$avg_log2FC > 0, ]
Symbol <- mapIds(get("org.Hs.eg.db"), keys = ors.markers$gene, keytype = "SYMBOL", column="ENTREZID")
ids <- bitr(ors.markers$gene,"SYMBOL","ENTREZID", "org.Hs.eg.db")
# 合并ENTREZID到sce.markers中
marker_data <- merge(ors.markers, ids, by.x="gene", by.y="SYMBOL")
gcSample <- split(marker_data$ENTREZID, marker_data$cluster)
gcSample
go <- compareCluster(gcSample, fun="enrichGO", OrgDb="org.Hs.eg.db", ont="BP", pvalueCutoff=0.05, qvalueCutoff=0.05)
res_go <- go@compareClusterResult
## 将富集结果中的 ENTREZID 重新转为 SYMBOL
for (i in 1:dim(res_go)[1]) {
  arr = unlist(strsplit(as.character(res_go[i,"geneID"]), split="/"))
  gene_names = paste(unique(names(Symbol[Symbol %in% arr])), collapse="/")
  res_go[i,"geneID"] = gene_names
}
head(res_go)
## 通路筛选
enrich_go_sbORS1 <- res_go %>% 
  group_by(Cluster) %>% 
  top_n(n = 15, wt = -pvalue) %>% 
  filter(Cluster %in% "sbORS1")
desired_GO_terms <- c("intermediate filament organization", "keratinization", "protein refolding", 
                      "cholesterol metabolic process", "skin development")
enrich_go_sbORS1 <- enrich_go_sbORS1 %>%
  filter(Description
         %in% desired_GO_terms)

enrich_go_sbORS2 <- res_go %>% 
  group_by(Cluster) %>% 
  top_n(n = 15, wt = -pvalue) %>% 
  filter(Cluster %in% "sbORS2")
desired_GO_terms <- c("axonogenesis", "establishment of cell polarity", "small GTPase-mediated signal transduction", 
                      "peptidyl-serine modification", "cell-cell signaling by wnt")
enrich_go_sbORS2 <- enrich_go_sbORS2 %>%
  filter(Description
         %in% desired_GO_terms)

enrich_go_bORS1 <- res_go %>% 
  group_by(Cluster) %>% 
  top_n(n = 15, wt = -pvalue) %>% 
  filter(Cluster %in% "bORS1")
desired_GO_terms <- c("response to ketone", "gland development", "transforming growth factor beta receptor superfamily signaling pathway", 
                      "regulation of BMP signaling pathway", "response to nutrient levels")
enrich_go_bORS1 <- enrich_go_bORS1 %>%
  filter(Description
         %in% desired_GO_terms)

enrich_go_bORS2 <- res_go %>% 
  group_by(Cluster) %>% 
  top_n(n = 15, wt = -pvalue) %>% 
  filter(Cluster %in% "bORS2")
desired_GO_terms <- c("cell-substrate adhesion", "Wnt signaling pathway", "small GTPase-mediated signal transduction", 
                      "cell-matrix adhesion", "axonogenesis")
enrich_go_bORS2 <- enrich_go_bORS2 %>%
  filter(Description
         %in% desired_GO_terms)

enrich_go_bORS3 <- res_go %>% 
  group_by(Cluster) %>% 
  top_n(n = 15, wt = -pvalue) %>% 
  filter(Cluster %in% "bORS3")
desired_GO_terms <- c("skin development", "epidermis development", "positive regulation of hydrolase activity", 
                      "keratinocyte differentiation", "regulation of response to wounding")
enrich_go_bORS3 <- enrich_go_bORS3 %>%
  filter(Description
         %in% desired_GO_terms)

enrich_go_bORS4 <- res_go %>% 
  group_by(Cluster) %>% 
  top_n(n = 15, wt = -pvalue) %>% 
  filter(Cluster %in% "bORS4")
desired_GO_terms <- c("stem cell differentiation", "establishment of cell polarity", "ribonucleoside diphosphate metabolic process", 
                      "regulation of GTPase activity", "ribonucleotide catabolic process")
enrich_go_bORS4 <- enrich_go_bORS4 %>%
  filter(Description
         %in% desired_GO_terms)


dt_sbORS1 <- enrich_go_sbORS1
dt_sbORS1 <- dt_sbORS1[order(dt_sbORS1$Cluster), ]
dt_sbORS1$Description <- factor(dt_sbORS1$Description, levels = dt_sbORS1$Description)
dt_sbORS1$geneID  <- sapply(strsplit(dt_sbORS1$geneID , "/"), function(x) paste(x[1:6], collapse = "/"))
dt_sbORS2 <- enrich_go_sbORS2
dt_sbORS2 <- dt_sbORS2[order(dt_sbORS2$Cluster), ]
dt_sbORS2$Description <- factor(dt_sbORS2$Description, levels = dt_sbORS2$Description)
dt_sbORS2$geneID  <- sapply(strsplit(dt_sbORS2$geneID , "/"), function(x) paste(x[1:6], collapse = "/"))
dt_bORS1 <- enrich_go_bORS1
dt_bORS1 <- dt_bORS1[order(dt_bORS1$Cluster), ]
dt_bORS1$Description <- factor(dt_bORS1$Description, levels = dt_bORS1$Description)
dt_bORS1$geneID  <- sapply(strsplit(dt_bORS1$geneID , "/"), function(x) paste(x[1:6], collapse = "/"))
dt_bORS2 <- enrich_go_bORS2
dt_bORS2 <- dt_bORS2[order(dt_bORS2$Cluster), ]
dt_bORS2$Description <- factor(dt_bORS2$Description, levels = dt_bORS2$Description)
dt_bORS2$geneID  <- sapply(strsplit(dt_bORS2$geneID , "/"), function(x) paste(x[1:6], collapse = "/"))
dt_bORS3 <- enrich_go_bORS3
dt_bORS3 <- dt_bORS3[order(dt_bORS3$Cluster), ]
dt_bORS3$Description <- factor(dt_bORS3$Description, levels = dt_bORS3$Description)
dt_bORS3$geneID  <- sapply(strsplit(dt_bORS3$geneID , "/"), function(x) paste(x[1:6], collapse = "/"))
dt_bORS4 <- enrich_go_bORS4
dt_bORS4 <- dt_bORS4[order(dt_bORS4$Cluster), ]
dt_bORS4$Description <- factor(dt_bORS4$Description, levels = dt_bORS4$Description)
dt_bORS4$geneID  <- sapply(strsplit(dt_bORS4$geneID , "/"), function(x) paste(x[1:6], collapse = "/"))

# Define the painting theme
mytheme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  axis.text.y = element_blank(), # 在自定义主题中去掉 y 轴通路标签:
  axis.ticks.length.y = unit(0,"cm"),
  plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5)
)
ggplot(dt_sbORS2, aes(x = -log10(qvalue), y = rev(Description), fill = Cluster))+
  geom_bar(stat = "identity", width = 0.5)+
  geom_text(aes(x=0.1,y=rev(Description),label = Description),size=6, hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(colour = 'black', size = 12),
        axis.line = element_line(colour = 'black', linewidth =0.5),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 12),
        legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_fill_manual(values = "#2682d1")+
  geom_text(data = dt_sbORS2,
            aes(x = 0.1, y = rev(Description), label = geneID, color = "#2682d1"),
            size = 4,
            fontface = 'italic', 
            hjust = 0,
            vjust = 4)+
  scale_color_manual(values ='#2682d1')+
  scale_y_discrete(expand = c(0.1,0))+
  labs(title = "GO Enrichment of BP",
       y= "sbORS2")


###   Rename   ###
Idents(ORS) <- factor(Idents(ORS), levels =  c("bORS1",
                                                     "bORS2",
                                                     "bORS3",
                                                     "bORS4",
                                                     "sbORS1",
                                                     "sbORS2" ) )
new.cluster.ids <- c("bORS-WI-SERPINF1","bORS-ECM-COL4A1", "bORS-SD-TNFRSF19","bORS-TD-PPARGC1A",
                     "sbORS-KR-KRT6B","sbORS-PL-PARD3")
names(new.cluster.ids) <- levels(ORS)
ORS <- RenameIdents(ORS, new.cluster.ids)
new.cluster.ids <- c("sbORS-KR-KRT6B","BG","bORS-SD-TNFRSF19", "sbORS-PL-PARD3","bORS-WI-SERPINF1",
                     "bTAC","CP","bORS-ECM-COL4A1","bORS-TD-PPARGC1A")
current.cluster.ids <- c("sbORS1","BG","bORS3","sbORS2","bORS1","bTAC","CP","bORS2","bORS4")
#替换分群序号为注释
ORS@meta.data$cell_type = plyr::mapvalues(x = ORS@meta.data[,"cell_type"], from = current.cluster.ids, to = new.cluster.ids)
new_order <- c("bORS-WI-SERPINF1","bORS-ECM-COL4A1", "bORS-SD-TNFRSF19","bORS-TD-PPARGC1A",
              "sbORS-KR-KRT6B","sbORS-PL-PARD3","CP","bTAC","BG")
# 重新定义 cell_type 的因子水平顺序
ORS@meta.data$cell_type <- factor(ORS@meta.data$cell_type, levels = new_order)

### ORS analysis

clusterCols <- c("#337ab7", "#4e99c7", "#a7cfe4","#d1e5f0",
                 "#21579d", "#2682d1")
names(clusterCols) <- c("bORS-WI-SERPINF1","bORS-ECM-COL4A1", "bORS-SD-TNFRSF19",  "bORS-TD-PPARGC1A",
                        "sbORS-KR-KRT6B","sbORS-PL-PARD3")

DimPlot(ORS, reduction = "umap", label = F, cols = clusterCols, pt.size = 0.3)+
  coord_fixed(ratio = 1)

ors.markers <- FindAllMarkers(ORS, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1.5) 
library(dplyr)
ors.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) 
top15 <- ors.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
bk <- c(seq(-3,-0.1),seq(0,3))
DoHeatmap(ORS, label = F, # 不加label
          features = as.character(unique(top15$gene)),   
          group.by = "cell_type",  
          assay = "RNA",  
          disp.min=-2,
          draw.lines = F,
          group.colors = c("#337ab7","#4e99c7","#a7cfe4","#d1e5f0","#21579d","#2682d1"))+ 
  scale_fill_gradientn(breaks=bk, colors = c("#93B5C1","white","#036862"))



###   Outer layer cells analysis   ###
## Rename
Idents(OL) <- factor(Idents(OL), levels =  c("bORS1",
                                              "bORS2",
                                              "bORS3",
                                              "bORS4",
                                              "sbORS1",
                                              "sbORS2",
                                              "CP",
                                              "bTAC",
                                              "BG") )
new.cluster.ids <- c("bORS-WI-SERPINF1","bORS-ECM-COL4A1", "bORS-SD-TNFRSF19","bORS-TD-PPARGC1A",
                     "sbORS-KR-KRT6B","sbORS-PL-PARD3", "CP","bTAC","BG")
names(new.cluster.ids) <- levels(OL)
OL <- RenameIdents(OL, new.cluster.ids)
new.cluster.ids <- c("sbORS-KR-KRT6B","BG","bORS-SD-TNFRSF19", "sbORS-PL-PARD3","bORS-WI-SERPINF1",
                     "bTAC","CP","bORS-ECM-COL4A1","bORS-TD-PPARGC1A")
current.cluster.ids <- c("sbORS1","BG","bORS3","sbORS2","bORS1","bTAC","CP","bORS2","bORS4")
#替换分群序号为注释
OL@meta.data$cell_type = plyr::mapvalues(x = OL@meta.data[,"cell_type"], from = current.cluster.ids, to = new.cluster.ids)
new_order <- c("bORS-WI-SERPINF1","bORS-ECM-COL4A1", "bORS-SD-TNFRSF19","bORS-TD-PPARGC1A",
               "sbORS-KR-KRT6B","sbORS-PL-PARD3","CP","bTAC","BG")
# 重新定义 cell_type 的因子水平顺序
OL@meta.data$cell_type <- factor(OL@meta.data$cell_type, levels = new_order)
clusterCols <- c("#337ab7", "#4e99c7", "#a7cfe4","#d1e5f0",
                 "#21579d", "#2682d1","#549da3","#96cb8f","#4dae47")
names(clusterCols) <- c("bORS-WI-SERPINF1","bORS-ECM-COL4A1", "bORS-SD-TNFRSF19",  "bORS-TD-PPARGC1A",
                        "sbORS-KR-KRT6B","sbORS-PL-PARD3", "CP","bTAC","BG")
DimPlot(OL, reduction = "umap", label = F, cols = clusterCols, pt.size = 0.3)+
  coord_fixed(ratio = 1)

# Pseudotime analysis
library(monocle3)
mat <- GetAssayData(OL, slot = "data")  
cellInfo <- OL@meta.data  
geneInfo <- data.frame(gene_short_name = rownames(mat), row.names = rownames(mat))  
cds <- new_cell_data_set(expression_data = mat,  
                         cell_metadata = cellInfo,  
                         gene_metadata = geneInfo)  
cds <- preprocess_cds(cds, num_dim = 50)
reducedDim(cds, "UMAP") <- Embeddings(OL, "umap")
cds = cluster_cells(cds, cluster_method = 'louvain')
cds = learn_graph(cds, use_partition=T, verbose=T, learn_graph_control=list(
  minimal_branch_len=10
)) #Trajectory learning

start = c("BG")
closest_vertex = cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex = as.matrix(closest_vertex[colnames(cds), ])
root_pr_nodes = igraph::V(principal_graph(cds)[["UMAP"]])$name
cell_ids <- rownames(colData(cds),"cell_type" == "BG")
flag = closest_vertex[cell_ids,]
flag = as.numeric(names(which.max(table(flag))))
root_pr_nodes = root_pr_nodes[flag]
cds = order_cells(cds, root_pr_nodes=root_pr_nodes)
##使用Seurat的UMAP信息，这样可以与Seurat对象的细胞分布保持一致
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(OL, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
clusterCols <- c("#337ab7", "#4e99c7", "#a7cfe4", "#d1e5f0",
                 "#21579d", "#2682d1", "#549da3", "#96cb8f", "#4dae47")
names(clusterCols) <- c("bORS-WI-SERPINF1", "bORS-ECM-COL4A1", "bORS-SD-TNFRSF19", 
                        "bORS-TD-PPARGC1A", "sbORS-KR-KRT6B", "sbORS-PL-PARD3", 
                        "CP", "bTAC", "BG")
plot_cells(cds, color_cells_by = "pseudotime",
           show_trajectory_graph=T) + plot_cells(cds,
                                                 color_cells_by = "cell_type",
                                                 label_cell_groups=FALSE,
                                                 label_leaves=FALSE,
                                                 label_branch_points=FALSE,
                                                 graph_label_size=1)+
  scale_color_manual(values = clusterCols) +
  scale_fill_manual(values = clusterCols)


#差异基因展示
Track_genes <- graph_test(cds,neighbor_graph="principal_graph", cores=6)
#按莫兰指数选择TOP基因
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
plot_genes_in_pseudotime(cds[Track_genes_sig,],color_cells_by="cell_type",min_expr=0.5, ncol= 2,cell_size=1.5) + 
  scale_color_manual(values = clusterCols)

genes <- c('SERPINF1','COL4A1','TNFRSF19','PPARGC1A','KRT6B','PARD3','KRT75','MKI67','CD34')
genes <- 'SERPINF1'
genes_cds <- cds[rowData(cds)$gene_short_name %in% genes,]
plot_genes_in_pseudotime(genes_cds,color_cells_by="cell_type",min_expr=0, ncol = 3, cell_size=1.5)+ 
  scale_color_manual(values = clusterCols)

###   Outer Layer Cell Ratio Analysis   ###
# Percent analysis
table(OL$orig.ident)#查看各组细胞数
prop.table(table(Idents(OL)))
table(Idents(OL), OL$orig.ident)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(OL), OL$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- data.frame(Cellratio)
library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]
###添加分组信息
sample <- c("F1","F2","F3","O1","O2","O3")
group <- c("Forehead","Forehead","Forehead","Occipital","Occipital","Occipital")
samples <- data.frame(sample, group)#创建数据框
rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']#R添加列
cellper$group <- samples[rownames(cellper),'group']#R添加列
###作图展示
pplist = list()
ol_groups = c("bORS-WI-SERPINF1", "bORS-ECM-COL4A1", "bORS-SD-TNFRSF19", 
               "bORS-TD-PPARGC1A", "sbORS-KR-KRT6B", "sbORS-PL-PARD3", 
               "CP", "bTAC", "BG")
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
for(group_ in ol_groups){
  cellper_  = cellper %>% dplyr::select(one_of(c('sample','group', group_)))#选择一组数据
  colnames(cellper_) = c('sample','group','percent')#对选择数据列命名
  cellper_$percent = as.numeric(cellper_$percent)#数值型数据
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))#上下分位数
  print(group_)
  print(cellper_$median)
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot作图
    geom_jitter(shape = 21,aes(fill=group),width = 0.25, size = 3) + 
    stat_summary(fun=mean, geom="point", color="grey60", size = 3) +
    theme_cowplot() +
    theme(axis.text = element_text(size = 7),axis.title = element_text(size = 7),legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),plot.title = element_text(size = 7,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
  ###组间t检验分析
  labely = max(cellper_$percent)
  compare_means(percent ~ group,  data = cellper_)
  my_comparisons <- list( c("Forehead", "Occipital") )
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons, method = "t.test", paired = T,size = 3)
  pplist[[group_]] = pp1
}
library(cowplot)
plot_grid(pplist[['bORS-WI-SERPINF1']],
          pplist[['bORS-ECM-COL4A1']],
          pplist[['bORS-SD-TNFRSF19']],
          pplist[['bORS-TD-PPARGC1A']],
          pplist[['sbORS-KR-KRT6B']],
          pplist[['sbORS-PL-PARD3']],
          pplist[['CP']],
          pplist[['bTAC']],
          pplist[['BG']]
)
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
ggplot(data = Cellratio, aes(x =Var2, y = Freq, fill =  Var1)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=c("#337ab7", "#4e99c7", "#a7cfe4", "#d1e5f0",
                             "#21579d", "#2682d1", "#549da3", "#96cb8f", "#4dae47")) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 






