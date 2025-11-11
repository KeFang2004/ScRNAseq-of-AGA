######################################
####     Fibroblasts Analysis     ####
######################################
FB <- subset(obj, idents = "FB")
FB <- ScaleData(FB, features = rownames(FB))
FB <- RunPCA(FB, features = VariableFeatures(object = FB),reduction.name = "pca")
FB <- RunHarmony(FB,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
FB <- RunUMAP(FB, reduction = "harmony", dims = 1:50,reduction.name = "umap")
DimPlot(FB, reduction = "umap",group.by = "group", label = T)
FB <- FindNeighbors(FB, reduction = "harmony", dims = 1:30)
FB_test <- FindClusters(FB, resolution = seq(0.1,1.9,0.2))
library(clustree)
clustree(FB_test@meta.data, prefix = "RNA_snn_res.")
rm(FB_test)
FB <- FindClusters(FB, resolution = 1.1)
FB <- identity(FB)
FB <- FB %>% RunTSNE(reduction = "harmony", dims = 1:30)
DimPlot(FB, reduction = "umap", group.by = "ident", split.by = "group", pt.size=1, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
FeaturePlot_scCustom(seurat_object = obj, features = "ACTA2",reduction = "umap",pt.size=0.1,
                     colors_use = c("#FFEFD5","#E6E6FA","#87CEFA","#6495ED","#4169E1","#0000CD","#000080"),min.cutoff =0.5)

names(FB@meta.data)
unique(FB$group)
Idents(FB) <- "group"
DEG_Findmarkers <- FindMarkers(FB,ident.1 = 'Forehead',ident.2 = 'Occipital', 
                               verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
DEG_Findmarkers <- DEG_Findmarkers %>%
  mutate(Type = if_else(p_val_adj > 0.05, "ns",
                        if_else(abs(avg_log2FC) < 1.5, "ns",
                                if_else(avg_log2FC >= 1.5, "up", "down")))) %>%
  arrange(desc(abs(avg_log2FC))) %>% rownames_to_column("Gene_Symbol")
table(DEG_Findmarkers$Type)
labeldata <- subset(DEG_Findmarkers, abs(avg_log2FC) > 1.5 & p_val_adj<0.05)
id <- order(-log10(labeldata$p_val_adj),decreasing = T)
labeldata <- labeldata[id,]
ggplot(DEG_Findmarkers,aes(avg_log2FC,-log10(p_val_adj)),max.overlaps = 100)+
  # 绘制散点
  geom_point(aes(color=Type),
             size=1.5)+
  # 绘制垂直线
  geom_vline(xintercept = c(-1.5,1.5), linetype = "dashed")+
  # 绘制水平线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  theme_bw()+
  # 设置点的颜色
  scale_color_manual(values = c('#5cb6ca','grey','#bb2e3d'))+
  # 设置标签注释
  geom_text_repel(data = labeldata[1:70,], 
                  aes(label = Gene_Symbol,color=Type),
                  max.overlaps = 70,
                  size=2)+
  # x-axis title
  xlab('Log2FC')+
  # y-axis title
  ylab('-Log10(adjPvalue)')
VlnPlot(FB, features = "FOSB",split.by = "group", cols = c('#bb2e3d','#5cb6ca'))

