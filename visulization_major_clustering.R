library(Seurat)
library(harmony)
library(ggplot2)
library(tidyverse)
library(stringr)
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

obj.markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 1.5) 
library(dplyr)
obj.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
write.csv(obj.markers, "primary_markers.csv")

my.colors <- colorRampPalette(c("lightblue", "white", "darkred"))(100)
pal<-viridis(n=15,option="D",direction=-1)
FeaturePlot_scCustom(seurat_object = obj, features = "KRT5",reduction = "umap",pt.size=0.1)
FeaturePlot(obj, features = "VCAN", reduction = "umap", pt.size = 0.1)

library(viridis)
library(scCustomize)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
FeaturePlot_scCustom(seurat_object = obj, features = gene,colors_use = c("#FFEFD5","#E6E6FA","#87CEFA","#6495ED","#4169E1","#0000CD","#000080"),
                     reduction = "umap",pt.size=0.01)+NoAxes()+
  theme(panel.border = element_rect(fill = NA,color = "black",
                                    size=1.5,linetype = "solid"))
Plot_Density_Joint_Only(seurat_object = obj, features = c("SFRP1", "MSX2"))
gene<-c("HLA-DQB2","CD3D","CD207","ACSBG1","AQP5","TAGLN","COL1A1","PECAM1","MLANA")
i=1
plots=list()
for (i in 1:length(gene)){
  plots[[i]]=FeaturePlot_scCustom(seurat_object = obj,
                                  colors_use = c("#FFEFD5","#E6E6FA","#87CEFA","#6495ED","#4169E1","#0000CD","#000080"),
                                  features = gene[i])+NoAxes()+
    theme(panel.border = element_rect(fill = NA,color = "black",
                                      size=1.5,linetype = "solid"))
}
library(patchwork)
p<-wrap_plots(plots, ncol = 3);p
p



Idents(obj) <- factor(Idents(obj), levels =  c("Mac","TC","LC",
                                                     "SC","GC","BC1","BC2","MC",
                                                     "SbG","SwG",
                                                     "UHF", "BG", "bORS", "sbORS", "CP", "TAC","HMC", "CO", "IRS", "FC", 
                                                     "ML","EC","FB","SMC"))
genes_to_check = c("HLA-DRB1", "IL32", "DOCK4", 
                   "IFI27", "CRABP2", "POSTN", "SAA1", "TK1",
                   "ACSBG1", "AQP5", 
                   "KRT17", "KRT15", "NEAT1", "SERPINA3", "CTSV", "MKI67", "MT1G", "SELENBP1","KRT28","KRT32",
                   "MLANA","PECAM1","COL1A1","TAGLN")

DotPlot(obj, features = genes_to_check, cluster.idents = F) +
  scale_color_gradientn(colours = c('#08306b', '#2170b5','#6baed6','#c6dbef','#fee0d2','#fc9271','#ef3b2c','#a50f15'))+
  theme(axis.text.x=element_text(angle = 90, hjust = 0.5,vjust=0.5))#+
#scale_y_discrete(limits = cell_order)  