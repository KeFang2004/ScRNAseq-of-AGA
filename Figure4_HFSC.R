library(Seurat)
library(harmony)
library(ggplot2)
library(tidyverse)
library(stringr)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(scCustomize)
library(monocle3)

FeaturePlot(BG, features = "CD34", reduction = "umap",pt.size = 0.3,order = T)
BG <- BG[, Idents(BG) %in% "BG"]

BG <- NormalizeData(BG, verbose = F) 
BG <- FindVariableFeatures(BG, selection.method = 'vst', nfeatures = 2000)
BG <- ScaleData(BG, vars.to.regress = "percent.mt")
BG <- RunPCA(BG, features = VariableFeatures(object = BG)) 
BG <- RunHarmony(BG,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
BG <- RunUMAP(BG, reduction = "harmony", dims = 1:30,reduction.name = "umap")
BG <- FindNeighbors(BG, dims = 1:30)
BG <- BG %>% RunTSNE(reduction = "harmony", dims = 1:30)
DimPlot(BG, reduction = "umap", label = TRUE, pt.size = 0.3)+
  coord_fixed(ratio = 1)
#BG_copy <- BG
#BG_copy <- FindClusters(BG_copy, resolution = seq(0.1,0.9,0.2))
#library(clustree)
#clustree(BG_copy@meta.data, prefix = "RNA_snn_res.")
BG <- FindClusters(BG, resolution = 0.3)
Idents(BG) <- 'seurat_clusters'
clusterCols <- c("#1ab17c", "#ea8bc3", "#93a8d8")
names(clusterCols) <- c("0","1","2")
DimPlot(BG, reduction = "umap", label = TRUE, cols = clusterCols, pt.size = 1)
 #+coord_fixed(ratio = 1)
DimPlot(BG, reduction = "umap", label = TRUE,pt.size = 1, split.by = "group")
FeaturePlot(BG, features = "CD200", reduction = "umap",pt.size = 0.3,order = T)

bg.markers <- FindAllMarkers(BG, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5) 
library(dplyr)
bg.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) 
top10 <- bg.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
FeaturePlot(BG, features = "CCDC3", reduction = "umap",pt.size = 0.3,order = T)
library(stringr)
library(viridis)
library(scCustomize)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggplot2)
FeaturePlot_scCustom(seurat_object = BG, features = "KRT7",reduction = "umap",
                     pt.size=1, colors_use = viridis_dark_high)

###############################################################
###   Violinplot showing differential expression of genes   ###
###############################################################
library(ggpubr)
library(ggimage)
library(ggplot2)
#设置比较-两两比较
my_comparisons <- list(c("Forehead", "Occipital"))
VlnPlot(BG, features = "WNT5A", group.by = "group")& 
  theme_bw()& 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(color = 'black',face = "bold", size = 12),
        axis.text.y = element_text(color = 'black', face = "bold"),
        axis.title.y = element_text(color = 'black', face = "bold", size = 15), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color="black",size = 1.2, linetype="solid"), 
        panel.spacing = unit(0.12, "cm"),
        plot.title = element_text(hjust = 0.5, face = "bold.italic"), 
        legend.position = 'none')&  
  stat_compare_means(method="t.test",hide.ns = F, comparisons = my_comparisons, 
                     label="p.signif",bracket.size=0.8, 
                     tip.length=0, size=6)&  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))&
  scale_fill_manual(values = c("#E5E4E2","#4A92A4"))





####################################
###   GO analysis for clusters   ###
####################################
# step1: Import the data #
options(stringsAsFactors = F)
library(dplyr)
table(bg.markers$cluster)
cluster0=bg.markers[bg.markers$cluster=='0',]$gene
cluster1=bg.markers[bg.markers$cluster=='1',]$gene
cluster2=bg.markers[bg.markers$cluster=='2',]$gene
total <- list(cluster0=cluster0,cluster1=cluster1,cluster2=cluster2)
# ID transfer
library(org.Hs.eg.db)
i=1
for (i in 1:3){
  #i=1
  ## 把SYMBOL改为ENTREZID
  total[[i]]=as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                        keys = total[[i]],
                                                        columns = 'ENTREZID',
                                                        keytype = 'SYMBOL')[,2]))}
lapply(total, head)
# step2: GO Enrichment #
library(clusterProfiler)
xx <- compareCluster(total, fun="enrichGO",
                     OrgDb = org.Hs.eg.db,
                     pvalueCutoff=0.95) # pvalueCutoff显著性应该为0.05
table(xx@compareClusterResult$Cluster)
head(as.data.frame(xx)) #Check the full result
xx <- simplify(xx,cutoff=0.4,by="p.adjust",select_fun=min)
p <- dotplot(xx) #气泡图
p + scale_color_gradient(low = "#d0ecd7", high = "#c1719c") +
  scale_fill_gradient(low = "#d0ecd7", high = "#c1719c")


df_go_diff <- as.data.frame(xx)
write.csv(df_go_diff, "GO_results_of_BG_clusters.csv")



#####################################
###   KEGG analysis for culsters  ###
#####################################
ids=bitr(bg.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
bg.markers=merge(bg.markers,ids,by.x='gene',by.y='SYMBOL')
gcSample=split(bg.markers$ENTREZID.x, bg.markers$cluster)
gcSample # entrez id , compareCluster 
KEGG_BG <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="hsa", 
                     pvalueCutoff= 1,
                     qvalueCutoff=1)
dotplot(KEGG_BG) 
# 生成 dotplot
p <- dotplot(KEGG_BG)

# 修改配色
# 使用 scale_color_gradient() 和 scale_fill_gradient() 修改颜色渐变
p + scale_color_gradient(low = "#cee9dd", high = "#3ba064") +
  scale_fill_gradient(low = "#cee9dd", high = "#3ba064")
p+ theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5))
KEGG_BG_df <- as.data.frame(KEGG_BG)
write.csv(KEGG_BG_df, "KEGG_results_of_BG_clusters.csv")

###############################
###   Difference Analysis   ###   
###############################
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(DESeq2)
library(AUCell)
library(reshape2)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)

bs = split(colnames(BG),BG$orig.ident)
a <- BG[["RNA"]]@features %>% rownames()
rownames(BG@assays$RNA@layers$counts) <- a
ct = do.call(
  cbind,lapply(names(bs), function(x){ 
    # x=names(bs)[[1]]
    kp =colnames(BG) %in% bs[[x]]
    rowSums( as.matrix(BG@assays$RNA@layers$counts[, kp]  ))
  })
)
colnames(ct) <- names(bs)
phe = unique(BG@meta.data[,c('orig.ident','group')])#样本&信息，自行修改
group_list = phe[match(names(bs),phe$orig.ident),'group']
exprSet = ct
exprSet=exprSet[apply(exprSet,1, function(x) sum(x>1) > 1),]
table(group_list)
#Step1 Build DEseq object
colData <- data.frame(row.names=colnames(exprSet),group_list=group_list)
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ group_list)
#Step2，进行差异表达分析
dds2 <- DESeq(dds)
table(group_list)
tmp <- results(dds2,contrast=c("group_list","Forehead","Occipital")) #分组自行修改
DEG_DESeq2 <- as.data.frame(tmp[order(tmp$padj),])
DEG_DESeq2 = na.omit(DEG_DESeq2)
#添加上下调信息
DEG_DESeq2 <- DEG_DESeq2 %>%
  mutate(Type = if_else(pvalue > 0.05, "ns",
                        if_else(abs(log2FoldChange) < 0.5, "ns",
                                if_else(log2FoldChange >= 0.5, "up", "down")))) %>%
  arrange(desc(abs(log2FoldChange))) %>% rownames_to_column("Gene_Symbol")
table(DEG_DESeq2$Type)

## 2)FindMarkers差异分析####
names(BG@meta.data)
unique(BG$group)
BG$celltype.group <- paste(BG$cell_type, BG$group, sep = "_")
Idents(BG) <- "group"
##选择某细胞亚群差异分析（以T cell为例）
#ident.1是要去比的组(Disease)，ident.2是被比较的组(Control)
DEG_Findmarkers <- FindMarkers(BG,ident.1 = 'Forehead',ident.2 = 'Occipital', 
                               verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
#添加上下调信息
DEG_Findmarkers <- DEG_Findmarkers %>%
  mutate(Type = if_else(p_val_adj > 0.05, "ns",
                        if_else(abs(avg_log2FC) < 0.5, "ns",
                                if_else(avg_log2FC >= 0.5, "up", "down")))) %>%
  arrange(desc(abs(avg_log2FC))) %>% rownames_to_column("Gene_Symbol")
table(DEG_Findmarkers$Type)


labeldata <- subset(DEG_Findmarkers, abs(avg_log2FC) > 0.5 & p_val_adj<0.05)
id <- order(-log10(labeldata$p_val_adj),decreasing = T)
labeldata <- labeldata[id,]

###   Volcano Plot   ###
ggplot(DEG_Findmarkers,aes(avg_log2FC,-log10(p_val_adj)))+
  # 绘制散点
  geom_point(aes(color=Type),
             size=1.5)+
  # 绘制垂直线
  geom_vline(xintercept = c(log2(1/1.4),log2(1.4)), linetype = "dashed")+
  # 绘制水平线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  theme_bw()+
  # 设置点的颜色
  scale_color_manual(values = c('#3ba064','#cee9dd','#FF857F'))+
  # 设置标签注释
  geom_text_repel(data = labeldata[1:40,], 
                  aes(label = Gene_Symbol,color=Type),
                  size=2.5)+
  # x-axis title
  xlab('Log2FC')+
  # y-axis title
  ylab('-Log10(adjPvalue)')

###############################
###   Enrichment Analysis   ###
###############################

# GSVA Analysis #
library(GSEABase)
library(GSVA)
geneSets <- getGmt('h.all.v2024.1.Hs.symbols.gmt')    ### Previously downloaded dataset
gsvaP <- ssgseaParam(
  exprData = exprSet,
  geneSets = geneSets,
  assay = NA_character_,
  annotation = NA_character_,
  minSize = 1,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE)
GSVA_hall <- gsva(gsvaP)
head(GSVA_hall)
GSVA_hall[,"O1"] <- GSVA_hall[,"O1"] + 0.2
# Visualization by limma #
library(limma)
# Import the group information
group <- factor(c(rep("Forehead", 3), rep("Occipital", 3)), levels = c('Forehead', 'Occipital'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(GSVA_hall)
design
# Forehead VS Occipital
compare <- makeContrasts(Forehead - Occipital, levels=design)
fit <- lmFit(GSVA_hall, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
head(Diff)
Diff <- rbind(subset(Diff,logFC>0)[1:10,], subset(Diff,logFC<0)[1:20,])
dat_plot <- data.frame(id = row.names(Diff),
                       t = Diff$t)
library(stringr)
dat_plot$id <- str_replace(dat_plot$id , "HALLMARK_","")
# 新增一列 根据t阈值分类
dat_plot$threshold = factor(ifelse(dat_plot$t  >-1, ifelse(dat_plot$t >= 1 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
# 排序
dat_plot <- dat_plot %>% arrange(t)
# Turn into the factors
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
library(ggplot2)
library(ggprism)
p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= '#36638a','NoSignifi'='#cccccc','Down'='#7bcd7b')) +
  geom_hline(yintercept = c(-1,1),color = 'white',linewidth = 0.5,lty='dashed') +
  xlab('') + 
  ylab('t value of GSVA score, Forehead versus Occipital') + #注意坐标轴旋转了
  guides(fill = F)+ # 不显示图例
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
# 小于-2的数量
low1 <- dat_plot %>% filter(t < -1) %>% nrow()
# 小于0总数量
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
# 小于2总数量
high0 <- dat_plot %>% filter(t < 1) %>% nrow()
# 总的柱子数量
high1 <- nrow(dat_plot)
p <- p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
                   hjust = 0,color = 'black') + # 小于-1的为黑色标签
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') +
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # 大于1的为黑色标签
p

#########################
###   GSEA analysis   ###
#########################

# 加载所需的R包
library(org.Hs.eg.db) # human的OrgDB
library(clusterProfiler)
# ID转化
gene_entrezid <- bitr(geneID = DEG_Findmarkers$Gene_Symbol, 
                      fromType = "SYMBOL", 
                      toType = "ENTREZID", # 转成ENTREZID
                      OrgDb = "org.Hs.eg.db"
)

head(gene_entrezid)
gene_entrezid$logFC <- DEG_Findmarkers$avg_log2FC[match(gene_entrezid$SYMBOL, DEG_Findmarkers$Gene_Symbol)]
genelist = gene_entrezid$logFC
names(genelist) = gene_entrezid$ENTREZID 
genelist=sort(genelist,decreasing = T) 
library(msigdbr)
m_t2g <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT_Legacy") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
gsea_res <- GSEA(genelist, 
                 TERM2GENE = m_t2g,
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 1,
                 pAdjustMethod = "BH"
)

library(enrichplot)
library(ggplot2)

ridgeplot(gsea_res,
          showCategory = 5,
          fill = "p.adjust", #填充色 "pvalue", "p.adjust", "qvalue" 
          core_enrichment = T,#是否只使用 core_enriched gene
          label_format = 30,
          orderBy = "p.adjust",
          decreasing = T
) +
  theme(axis.text.y = element_text(size=8))

ids <- gsea_res@result$ID[1:5]
gseadist(gsea_res,
         IDs = ids,
         type = "density" # boxplot
) +
  theme(legend.direction = "vertical")

ridgeplot(gsea_res,
          showCategory = 5,
          fill = "p.adjust", # 填充色 "pvalue", "p.adjust", "qvalue"
          core_enrichment = T, # 是否只使用 core_enriched gene
          label_format = 30,
          orderBy = "p.adjust",
          decreasing = T
) +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Adjusted\nP-value") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major.y = element_line(color = "gray80"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  ) +
  labs(
    title = "GSEA Analysis Ridge Plot",
    x = "Enrichment Score",
    y = "Pathway",
    fill = "Adjusted P-value"
  )


GSEA_df <- as.data.frame(gsea_res)
write.csv(GSEA_df, "Figure4_HFSC_TF_GSEA.csv")
library(enrichplot) # 富集结果可视化
# 特定通路作图——单个通路
gseaplot2(gsea_res, "ATF3_Q6", color = "steelblue", pvalue_table = T)


#######网络图########
## 将富集结果中的 ENTREZID 重新转为 SYMBOL
Symbol <- mapIds(get("org.Hs.eg.db"), keys = DEG_Findmarkers$Gene_Symbol, keytype = "SYMBOL", column="ENTREZID")
for (i in 1:dim(GSEA_df)[1]) {
  arr = unlist(strsplit(as.character(GSEA_df[i,"core_enrichment"]), split="/"))
  gene_names = paste(unique(names(Symbol[Symbol %in% arr])), collapse="/")
  GSEA_df[i,"core_enrichment"] = gene_names
}
gsea_res1 <- setReadable(gsea_res, 'org.Hs.eg.db', 'ENTREZID')
p <- cnetplot(gsea_res1, showCategory =5, categorySize = "pvalue", colorEdge = T,foldChange=genelist, color_category ='#FF857F')
# 运行cnetplot，并设置颜色参数
p <- cnetplot(gsea_res1, 
              foldChange = genelist,
              layout = 'drl',
              color.params = list(foldChange = genelist,
                                  edge = TRUE),
              categorySize = "pvalue",
              node_label = "all")
p + scale_color_gradientn(colours = c("#6d3e99", "#f4e2ee"))





############################
###   Heatmap Painting   ###
############################
library(ComplexHeatmap)
library(pheatmap)
Expmatrix <- GetAssayData(BG, layer = "scale.data")
IL2_genes <- c(
  "ABCB1", "ADAM19", "AGER", "AHCY", "AHNAK", "AHR", "AKAP2", "ALCAM", "AMACR", 
  "ANXA4", "APLP1", "ARL4A", "BATF", "BATF3", "BCL2", "BCL2L1", "BHLHE40", "BMP2", 
  "BMPR2", "CA2", "CAPG", "CAPN3", "CASP3", "CDC164", "CCND2", "CCND3", "CCNE1", 
  "CCR4", "CD44", "CD48", "CD79B", "CD81", "CD83", "CD86", "CDC42SE2", "CDC6", 
  "CDCP1", "CDKN1C", "CISH", "CKAP4", "COCH", "COL6A1", "CSF1", "CSF2", "CST7", 
  "CTLA4", "CTSZ", "CXCL10", "CYFIP1", "DCPS", "DENND5A", "DHRS3", "ECM1", "EMP1", 
  "ENO3", "ENPP1", "EOMES", "ETV4", "F2RL2", "FAH", "FAM126B", "FGL2", "FLT3LG", 
  "FURIN", "GABARAPL1", "GADD45B", "GALM", "GATA1", "GBP4", "GLIPR2", "GPR65", 
  "GPR83", "GPX4", "GSTO1", "GUCY1B3", "HIPK2", "HK2", "HOPX", "HUWE1", "ICOS", 
  "IFITM3", "IFNGR1", "IGF1R", "IGF2R", "IKZF2", "IKZF4", "IL10", "IL10RA", "IL13", 
  "IL18R1", "IL1R2", "IL1RL1", "IL2RA", "IL2RB", "IL3RA", "IL4R", "IRF4", "IRF6", 
  "IRF8", "ITGA6", "ITGAE", "ITGAV", "ITIH5", "KLF6", "LCLAT1", "LIF", "LRIG1", 
  "LRRC8C", "LTB", "MAFF", "MAP3K8", "MAP6", "MAPKAPK2", "METTL20", "MUC1", "MXD1", 
  "MYC", "MYO1C", "MYO1E", "N6AMT2", "NCOA3", "NCS1", "NDRG1", "NFIL3", "NFKBIZ", 
  "NOP2", "NRP1", "NT5E", "ODC1", "P2RX4", "P4HA1", "PDCD2L", "PENK", "PHLDA1", 
  "PHTF2", "PIM1", "PLAGL1", "PLEC", "PLIN2", "PLSCR1", "PNP", "POU2F1", "PPAP2A", 
  "PRAF2", "PRKCH", "PRNP", "PTCH1", "PTGER2", "PTH1R", "PTRH2", "PUS1", "RABGAP1L", 
  "RGS16", "RHOB", "RHOH", "RNH1", "RORA", "RRAGD", "S100A1", "SCN9A", "SELL", 
  "SELP", "SERPINB6", "SERPINC1", "SH3BGRL2", "SHE", "SLC1A5", "SLC29A2", "SLC2A3", 
  "SLC39A8", "SMPDL3A", "SNX14", "SNX9", "SOCS1", "SOCS2", "SPP1", "SPRED2", "SPRY4", 
  "ST3GAL4", "SWAP70", "SYNGR2", "SYT11", "TGM2", "TIAM1", "TLR7", "TNFRSF18", 
  "TNFRSF1B", "TNFRSF21", "TNFRSF4", "TNFRSF8", "TNFRSF9", "TNFSF10", "TNFSF11", 
  "TRAF1", "TTC39B", "TWSG1", "UCK2", "UMPS", "WLS", "XBP1")
WNT_genes <-c(
  "ADAM17", "AXIN1", "AXIN2", "CCND2", "CSNK1E", "CTNNB1", "CUL1", "DKK1", "DKK4", 
  "DLL1", "DVL2", "FZD1", "FZD8", "GNAI1", "HDAC11", "HDAC2", "HDAC5", "HCEY1", 
  "HEY2", "JAG1", "JAG2", "KAT2A", "LEF1", "AML1", "MYC", "NCOR2", "NCSTN", 
  "NKD1", "NOTCH1", "NOTCH4", "NUMB", "PPARD", "PSEN2", "PTCH1", "RBPJ", "SKP2", 
  "TCF7", "TP53", "WNT1", "WNT5", "WNT6")
TNF_genes <- c(
  "ABCA1", "ACKR3", "AREG", "ATF3", "ATP2B1", "B4GALT1", "B4GALT5", "BCL2A1", "BCL3", "BCL6", 
  "BHLHE40", "BIRC2", "BIRC3", "BMP2", "BTG1", "BTG2", "BTG3", "CCL2", "CCL20", "CCL4", "CCL5", 
  "CCN1", "CCND1", "CCNL1", "CCRL2", "CD44", "CD69", "CD80", "CD83", "CDKN1A", "CEBPB", "CEBPD", 
  "CFLAR", "CLCF1", "CSF1", "CSF2", "CXCL1", "CXCL10", "CXCL11", "CXCL2", "CXCL3", "CXCL6", "DENND5A", 
  "DNAJB4", "DRAM1", "DUSP1", "DUSP2", "DUSP4", "DUSP5", "EDN1", "EFNA1", "EGR1", "EGR2", "EGR3", 
  "EHD1", "ETS2", "F2RL", "F3", "FJX1", "FOSB", "FOSL2", "FUT4", "G0S2", "GADD45A", "GADD45B", "GEM", 
  "GFPT2", "GPR183", "HBEGF", "HES1", "ICAM1", "ICOSLG", "IDER2", "IER3", "IER5", "IFIH1", "IFIT2", 
  "IFNGR2", "IL12B", "IL15RA", "IL18R1", "IL1B", "IL2A6", "IL6ST", "IL7R", "INHBA", "IRF1", "IRS2", 
  "JAG1", "JUNB", "KDM6B", "KLF10", "KLF2", "KLF4", "KLF6", "KLF9", "KYNU", "LAMB3", "LDLDR", "LIF", 
  "LITAF", "MAFF", "MAP2K3", "MAP3K8", "MARCKS", "MCL1", "MSC", "MXD1", "MYC", "NAMPT", "NFAT5", "NFE2L2", 
  "NFIL3", "NFKB1", "NFKB2", "NFKBIA", "NFKBIE", "NINJ1", "NR4A1", "NR4A2", "NR4A3", "OLR1", "PANX1", 
  "PDE4B", "PDLIM5", "PER1", "PFKFB3", "PHLDA1", "PHLDA2", "PLAU", "PLEK", "PLK2", "PLPP3", "PMEPA1", 
  "PNRC1", "PPR15A", "PTGER4", "PTGS2", "PTX3", "RCAN1", "RELA", "RHOB", "RIGI", "RIPK2", "RNF19B", 
  "S100A", "SERPIN2", "SERPINB8", "SINE1", "SGK1", "S1K1", "SLC16A6", "SLC2A3", "SLC2A6", "SMAD3", 
  "SNN", "SOCS3", "SOD2", "SPHK1", "SPB1", "SQSTM1", "STAT5A", "TANK")
DNA_gene <- c(
  "AAAS", "ADA", "ADCY6", "ADRM1", "AGO4", "AK1", "AK3", "ALYREF", "ART", "ARL6IP1", 
  "BCAM", "BCAP31", "BOLA2", "BRF2", "CANT1", "CCNO", "CDA", "CTN2", "CLP1", "CMPK2", 
  "COX17", "CSTF3", "DAD1", "DCTN4", "DDB1", "DGB2", "DGUOK", "DUT", "EDF1", 
  "EIFB", "ELLA", "ERCC1", "ERC2", "ERC3", "ERC4", "ERC5", "ERC8", "FEN1", "GMPR2", 
  "GPX4", "GSDME", "GTF2A2", "GTF2B", "GTF2F1", "GTF2H1", "GTF2H3", "GTF2H5", 
  "GTF3C5", "GUK1", "HCLS1", "HPR", "IMPDH2", "ITPA", "LIG1", "MPC2", "MPG", 
  "MRPL0", "NCBP2", "NELFB", "NELFCD", "NFX1", "NME1", "NME3", "NME4", "NPR2", 
  "NT5C", "NT5C3A", "NUDT2", "NUDT9", "PCNA", "PDE4B", "PDE6G", "PNPNA", "POLA2", 
  "POLB", "POLD1", "POLE4", "POLL", "POLR1C", "POLR1D", "POLR1H", "POLR2A", "POLR2C", 
  "POLR2D", "POLR2E", "POLR2F", "POLR2G", "POLR2H", "POLR2I", "POLR2J", "POLR2K", "POLR3C", 
  "POLR3GL", "POM121", "PMRIM1", "RAD51", "RAD52", "RAE1", "RALA", "RBX1", "REV3L", "RFC2", 
  "RFC3", "RFC4", "RFC5", "RNMT", "RPA2", "RPA3", "RM2B", "SAC3D1", "SDCBP", "SEC6A1", 
  "SF3A3", "SMA", "SNAPC4", "SNAPC5", "SRSF6", "SRP1", "STX3", "SUPT4H1", "SUPT5H", "SURF1", 
  "TAF10", "TAF12", "TAF13", "TAF1C", "TAF6", "TAF9", "TARBP2", "TK2", "TMED2", "TP53", 
  "TSG101", "TYMS", "UMPS", "UPF3B", "USP11", "VPS28", "VPS37B", "XPC", "ZNF707", "ZWINT"
)
cluster_info <- sort(BG$orig.ident)
# 获取Expmatrix中实际存在的基因名
existing_genes <- rownames(Expmatrix)
# 筛选出在Expmatrix中存在的IL2基因
valid_genes <- DNA_gene[DNA_gene %in% existing_genes]
Expmatrix<- as.matrix(Expmatrix[valid_genes,names(cluster_info)])
col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
col_fun(seq(-3, 3))
Heatmap(Expmatrix, column_split = cluster_info)
Heatmap(Expmatrix, column_split = cluster_info, col = col_fun)
row_mean = apply(Expmatrix,1,mean)


gene_cell_exp <- AverageExpression(BG,
                                   features = valid_genes,
                                   group.by = 'orig.ident',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
#顶部细胞类型注释
df <- data.frame(colnames(gene_cell_exp))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,#细胞名/cluster
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('F1'="#E5E4E2",
                                                  'F2'="#E5E4E2",
                                                  'F3'="#E5E4E2",
                                                  'O1'="#4A92A4",
                                                  'O2'="#4A92A4",
                                                  'O3'="#4A92A4")))#颜色设置
#数据标准化缩放一下
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
Heatmap(marker_exp,
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = T,
        column_title = NULL,
        heatmap_legend_param = list(
          title=' '),
        #col = colorRampPalette(c("#d0ecd7","white","#c1719c"))(100),
        border = 'black',
        rect_gp = gpar(col = "black", lwd = 0.1),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        top_annotation = top_anno)

###########################
###   SCENIC Analysis   ###
###########################

## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
#BiocManager::install(c("zoo", "mixtools", "rbokeh"),ask = F,update = F) 
# For various visualizations and perform t-SNEs:
#BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"),ask = F,update = F)
# To support paralell execution (not available in Windows):
#BiocManager::install(c("doMC", "doRNG"),ask = F,update = F)
#devtools::install_github("aertslab/SCENIC") 
install.packages("arrow") #正常情况下由于网络原因会安装不完全
library(arrow)   #导致加载包有警告，按提示输入一下内容
install_arrow(verbose=TRUE)
library(SCENIC)
exprMat = GetAssayData(BG, slot="counts")
cellInfo <-  BG@meta.data[,c(4,2,3)]
colnames(cellInfo)=c('CellType', 'nGene' ,'nUMI')
head(cellInfo)
table(cellInfo$CellType)
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
data(defaultDbNames)
dbs <- defaultDbNames[['hgnc']]
dbs[[1]] <- "hg19-500bp-upstream-7species.mc9nr.feather"
dbs[[2]] <-"hg19-tss-centered-10kb-7species.mc9nr.feather"
scenicOptions <- initializeScenic(org="hgnc", dbDir="~/AGA_ENDO_RE", dbs = dbs, nCores=1) 

