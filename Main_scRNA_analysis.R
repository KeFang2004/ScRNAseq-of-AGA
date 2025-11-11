library(Seurat)
library(harmony)
library(DoubletFinder)
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
library(simplifyEnrichment)
library(plot1cell)
library(patchwork)
library(Matrix)
library(biomaRt)

Forehead = "Forehead"
Occipital = "Occipital"
F1 <- CreateSeuratObject(counts = Read10X('F-1'), 
                         project = "Forehead", 
                         min.cells = 5, 
                         min.features = 300)
F1$orig.ident <- "F1"
F1$group <- "Forehead"

F2 <- CreateSeuratObject(counts = Read10X('F-2'), 
                         project = "Forehead", 
                         min.cells = 5, 
                         min.features = 300)
F2$orig.ident <- "F2"
F2$group <- "Forehead"

F3 <- CreateSeuratObject(counts = Read10X('F-3'), 
                         project = "Forehead", 
                         min.cells = 5, 
                         min.features = 300)
F3$orig.ident <- "F3"
F3$group <- "Forehead"

O1 <- CreateSeuratObject(counts = Read10X('OC-1'), 
                         project = "Occipital", 
                         min.cells = 5, 
                         min.features = 300)
O1$orig.ident <- "O1"
O1$group <- "Occipital"

O2 <- CreateSeuratObject(counts = Read10X('OC-2'), 
                         project = "Occipital", 
                         min.cells = 5, 
                         min.features = 300)
O2$orig.ident <- "O2"
O2$group <- "Occipital"

O3 <- CreateSeuratObject(counts = Read10X('OC-3'), 
                         project = "Occipital", 
                         min.cells = 5, 
                         min.features = 300)
O3$orig.ident <- "O3"
O3$group <- "Occipital"

# Put the Seurat object into a list
data <- list(F1, F2, F3, O1, O2, O3)

# Combine the Seurat list and keep the sample names
obj <- merge(x = data[[1]], y = data[-1], 
              project = "AGA")

# Check the sample names and group information
head(obj$orig.ident)
head(obj$group)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), alpha = 0.4,ncol = 3) #Quality control visualization

obj <- subset(obj, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
obj <- NormalizeData(obj, verbose=F)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(obj), 10)                                         #Show top 10 variable genes

#plot1 <- VariableFeaturePlot(obj)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1 + plot2

obj <- ScaleData(obj, features = rownames(obj))
obj <- RunPCA(obj, features = VariableFeatures(object = obj),reduction.name = "pca")

obj <- RunHarmony(obj,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30,reduction.name = "umap")
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)

### Find the best resolution using clustree package
#obj_test <- FindClusters(obj, resolution = seq(0.1,1.3,0.2))
#library(clustree)
#clustree(obj_test@meta.data, prefix = "RNA_snn_res.")
# The best resolution we found was 0.5
obj <- FindClusters(obj, resolution = 0.5)
obj <- identity(obj)
# Now there are 24 clusters

### Remove doublets
sweep.res.list <- paramSweep(obj, PCs = 1:20, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) 
sweep.stats[order(sweep.stats$BCreal),]
bcmvn <- find.pK(sweep.stats) 
pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)]) 
homotypic.prop <- modelHomotypic(scRNA_harmony$seurat_clusters) 
nExp_poi <- round(0.075 *nrow(scRNA_harmony@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) 
# Identify doublets
scRNA_harmony <- doubletFinder_v3(scRNA_harmony, PCs = 1:20, pN = 0.25, pK = pK_bcmvn,
                                  nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
obj <- subset(scRNA_harmony, subset = DF.classifications_0.25_33_9717== "Singlet")

# TSNE
obj <- obj %>% 
  RunTSNE(reduction = "harmony", dims = 1:30)
# See UMAP reduction result 
DimPlot(obj, reduction = "umap", group.by = "ident", split.by = "orig.ident", pt.size=0.5)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
# See cluster result
DimPlot(obj, reduction = "umap", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)

# Find cell type marker genes
obj <- JoinLayers(obj)
cluster.markers <- FindMarkers(obj, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
obj.markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1.5) 
library(dplyr)
obj.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) 
top10 <- obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
VlnPlot(obj, features = top10$gene[1:20],pt.size=0)

# Try to use previous verified marker genes to annotate the clusters
# The markers we used can be found in the references of the manuscript, which were verified by previous researches
genes_to_check = list(
  BC = c("KRT5", "KRT14", "KRT15", "COL17A1", "CDH3"),
  MC = c("STMN1", "HMGB2", "MKI67"),
  SC = c("KRT1", "KRT10", "MT4", "KRT77"),
  GC = c("CRCT1", "FLG2"),
  UHF = c("KRT17", "CCL27A", "SOX9", "KRT79"),
  HFB = c("CTGF", "CXCL14", "POSTN"),
  GL = c("SOX4", "PTN", "FABP5"),
  SG = c("SCD", "MGST1", "CERS4"),
  EC = c("FABP4", "PECAM1"),
  SM = c("TAGLN","ACTA2"),
  FB = c("COL1A1", "COL3A1", "FN1"),
  Mac = c("CD74","C1QA"),
  LC = c("CD207", "CD52", "CD86"),
  TC = c("CTLA2A", "CD3G", "CD7")
)
genes_to_check = c("KRT14", "MKI67", "KRT1", "FLG2", "KRT79", "LHX2", "CD34", "CDH3", "LEF1", "SCD", "AQP5", "CD68", "CD163", "CD207", "CD3G", "PECAM1", "ACTA2", "PDGFRA")
DotPlot(obj,group.by = 'seurat_clusters', features = genes_to_check,cluster.idents = T) +
  scale_color_gradientn(colours = c('#FFCC33', '#66CC66','#336699','#330066'))+
  theme(axis.text.x=element_text(angle = 90, hjust = 0.5,vjust=0.5))#+
  #scale_y_discrete(limits = cell_order)  
FeaturePlot(obj, features = "KRT17", reduction = "umap", pt.size = 0.1)
VlnPlot(obj,features = "CST6")
library(Nebulosa)
library(viridisLite)
library(viridis)
library(scCustomize)
Plot_Density_Joint_Only(obj, reduction = "umap",
                   features = c("KRT5", "BARX2"))
# GO enrichment analysis to identify the function of different clusters
library(clusterProfiler)
top100 <- obj.markers %>% group_by(cluster) %>% top_n(100, avg_log2FC)
de_genes <- subset(top100, p_val_adj<0.05)
length(de_genes$gene) # 1000
head(de_genes)
# Transfer gene ID
entrez_genes <- bitr(de_genes$gene, fromType="SYMBOL", 
                     toType="ENTREZID", 
                     OrgDb="org.Hs.eg.db") 
head(entrez_genes)
# Match the previous gene symbol with the entrez id
mix <- merge(de_genes,entrez_genes,by.x="gene",by.y="SYMBOL")
dim(mix) 
head(mix)
de_gene_clusters <- data.frame(ENTREZID=mix$ENTREZID,
                               cluster=mix$cluster)
table(de_gene_clusters$cluster)
formula_res <- compareCluster(
  ENTREZID~cluster, 
  data=de_gene_clusters, 
  fun="enrichGO", 
  OrgDb="org.Hs.eg.db",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  =0.05,
  minGSSize = 5, 
  maxGSSize = 1000) 
GO = as.data.frame(formula_res)
write.csv(GO, "GO_analyasis_for_seurat_cluster.csv")

cluster2celltype <- c("0"="SC",
                      "1"="BC1", 
                      "2"="sbORS", 
                      "3"= "GC", 
                      "4"= "BG", 
                      "5"= "MC",
                      "6"= "CO", 
                      "7"= "UHF", 
                      "8"= "GC",
                      "9"= "BC2",
                      "10"= "bORS",
                      "11"= "HMC",
                      "12"= "bORS",
                      "13"= "ML",
                      "14"= "TAC",
                      "15"= "IRS",
                      "16"= "IC",
                      "17"= "EC",
                      "18"= "FB",
                      "19"= "IC",
                      "20"= "CP",
                      "21"= "SwG",
                      "22"= "FC",
                      "23"= "SMC",
                      "24"= "SbG")
new.cluster.ids <- c("SC", "BC1", "sbORS", "GC", "BG", "MC", "CO", "UHF", "GC", "BC2", "bORS","HMC",
                     "bORS","ML","TAC","IRS","IC","EC","FB","IC","CP","SwG","FC","SMC","SbG")
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj, new.cluster.ids)
current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)
#Change the number order into the annotations
obj@meta.data$cell_type = plyr::mapvalues(x = obj@meta.data[,"seurat_clusters"], from = current.cluster.ids, to = new.cluster.ids)
obj <- obj[, order(match(obj@meta.data$merge_celltype, cell_order))]

#clusterCols <- c("#85C4C9", "#386DAA", "#B0F2BC",
#                 "#FBE6C5", "#78A02D", "#96D2A4", 
#                 "#3099B1", "#FF7AA3", "#B099C9", "#A4A494",
#                 "#9D1A42", 
#                 "#B185DC", "#FDE086", "#F46C44")
#names(clusterCols) <- c("SC", "BC", "UHF", "SbG", "HFB", "GL", "MC", "Mac", "TC", "LC","ML", "FB", "EC","SwG")
#DimPlot(obj, reduction = "umap", cols = clusterCols, label = TRUE, pt.size = 0.1)
#DimPlot(obj, reduction = "umap", split.by = "group", label = TRUE, pt.size = 0.1)

library(plot1cell)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(biomaRt)
complex_dotplot_single(seu_obj = obj, feature = "KRT17",groups = "orig.ident")


# As we know, immune cells in the micro-environment of hair follicles can defenitely be divided into specific clusters, such as macrophages and T cells
# Divide immune cells from IC 
sub_cells <- subset(obj, cell_type == "IC")
# Determine the K-nearest neighbor graph
sub_cells <- FindNeighbors(sub_cells, reduction="pca", dims = 1:20, verbose = F)
# Determine the clusters for various resolutions
sub_cells <- FindClusters(sub_cells, resolution = 0.1, verbose = F)
# t-SNE and Clustering
sub_cells <- RunTSNE(sub_cells, reduction = "pca", dims = 1:20, verbose = F) # check_duplicates = FALSE, Remove duplicates before running TSNE
sub_cells <- RunUMAP(sub_cells, reduction = "pca", dims = 1:20, verbose = F)
DimPlot(sub_cells, reduction = "umap", split.by = "group", label = TRUE, pt.size = 0.1)
genes_to_check = c("CD68", "CD207","CD3G")
DotPlot(sub_cells,group.by = 'seurat_clusters', features = genes_to_check,cluster.idents = T) +
  scale_color_gradientn(colours = c('#FFCC33', '#66CC66','#336699','#330066'))+
  theme(axis.text.x=element_text(angle = 90, hjust = 0.5,vjust=0.5))
cluster2celltype <- c("0"="Mac",
                      "1"="TC", 
                      "2"="LC"
                   )
new.cluster.ids <- c("Mac", "TC", "LC")
names(new.cluster.ids) <- levels(sub_cells)
sub_cells <- RenameIdents(sub_cells, new.cluster.ids)
sub_cells[['cell_type']] = unname(cluster2celltype[sub_cells@meta.data$seurat_clusters])
sub_cells@meta.data$seurat_clusters <-  sub_cells@meta.data$celltype 
sub_cells$new_celltype <- sub_cells$cell_type
## Mapping to old clusters
Idents(obj, cells = colnames(sub_cells)) <- Idents(sub_cells)
DimPlot(obj,label = TRUE)


clusterCols <- c("#96cb8f","#8bc96d","#4dae47","#5c9e43","#a3e0bb","#BCECD1",
                 "#f9b769","#f9a341","#d480a6a8","#c6b598",
                 "#a4cde1","#67a4cc","#277fb8","#549da3",
                 "#f6f28f","#d4a55b","#b05a28",
                 "#F4B37A","#af8093c4","#EFA48F","#f38989","#9886A7","#815e99","#EDB78E")
names(clusterCols) <- c("sbORS", "UHF", "BG", "TAC", "CP", "bORS", 
                        "SMC", "FB", "EC", "ML",
                        "HMC", "IRS", "CO","FC",
                        "TC","Mac","LC",
                        "BC1","SwG","BC2","GC","MC","SbG","SC")
DimPlot(obj, reduction = "umap",cols = clusterCols,label = TRUE, pt.size = 0.1)


###   Percent analysis
table(obj$orig.ident)#Check cell numbers in two groups
prop.table(table(Idents(obj)))
table(Idents(obj), obj$orig.ident) #Cell numbers of different cell types in each group
Cellratio <- prop.table(table(Idents(obj), obj$orig.ident), margin = 2) # Calclulate the cell ratio
Cellratio
Cellratio <- data.frame(Cellratio)
library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq") 
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]
# Add the grouping information
sample <- c("F1","F2","F3","O1","O2","O3")
group <- c("Forehead","Forehead","Forehead","Occipital","Occipital","Occipital")
samples <- data.frame(sample, group) #Create a data frame
rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']# Add columns
cellper$group <- samples[rownames(cellper),'group']# Add columns
# Draw
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
pplist = list()
sce_groups = c("sbORS", "UHF", "BG", "TAC", "CP", "bORS", 
               "SMC", "FB", "EC", "ML",
               "HMC", "IRS", "CO","FC",
               "TC","Mac","LC",
               "BC1","SwG","BC2","GC","MC","SbG","SC")

for(group_ in sce_groups){
  cellper_  = cellper %>% dplyr::select(one_of(c('sample','group', group_)))# Select a group of the data
  colnames(cellper_) = c('sample','group','percent')# Name
  cellper_$percent = as.numeric(cellper_$percent) # Change into the numeric data
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent)) 
  print(group_)
  print(cellper_$median)
  
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + 
    geom_jitter(shape = 21,aes(fill=group),width = 0.25, size = 3) + 
    stat_summary(fun=mean, geom="point", color="grey60", size = 3) +
    theme_cowplot() +
    theme(axis.text = element_text(size = 7),axis.title = element_text(size = 7),legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),plot.title = element_text(size = 7,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
  
  ### t-test between the two groups
  labely = max(cellper_$percent)
  compare_means(percent ~ group,  data = cellper_)
  my_comparisons <- list( c("Forehead", "Occipital") )
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons, method = "t.test", size = 3)
  pplist[[group_]] = pp1
}


library(cowplot)
sce_groups = c("sbORS", "UHF", "BG", "TAC", "CP", "bORS", 
               "SMC", "FB", "EC", "ML",
               "HMC", "IRS", "CO","FC",
               "TC","Mac","LC",
               "BC1","SwG","BC2","GC","MC","SbG","SC")
plot_grid(pplist[['sbORS']],
          pplist[['UHF']],
          pplist[['BG']],
          pplist[['TAC']],
          pplist[['CP']],
          pplist[['bORS']],
          pplist[['SMC']],
          pplist[['FB']],
          pplist[['EC']],
          pplist[['ML']],
          pplist[['HMC']],
          pplist[['IRS']],
          pplist[['CO']],
          pplist[['FC']],
          pplist[['TC']],
          pplist[['Mac']],
          pplist[['LC']],
          pplist[['BC1']],
          pplist[['SwG']],
          pplist[['BC2']],
          pplist[['GC']],pplist[['MC']],pplist[['SbG']],pplist[['SC']]
          )
Cellratio <- as.data.frame(Cellratio)
write.csv(Cellratio,"Cellratio_ALL.csv")
colourCount = length(unique(Cellratio$Var1))
ggplot(data = Cellratio, aes(x =Var2, y = Freq, fill =  Var1)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=c("#96cb8f","#8bc96d","#4dae47","#5c9e43","#a3e0bb","#BCECD1",
                             "#f9b769","#f9a341","#d480a6a8","#c6b598",
                             "#a4cde1","#67a4cc","#277fb8","#549da3",
                             "#f6f28f","#d4a55b","#b05a28",
                             "#F4B37A","#af8093c4","#EFA48F","#f38989","#9886A7","#815e99","#EDB78E")) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  #### Move the location of the axis
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 

###   Check gene expression level 
library(ggpubr)
library(ggimage)
library(ggplot2)
my_comparisons <- list(c("Forehead", "Occipital"))
VlnPlot(obj, features = "OVOL1", group.by = "group", idents = "HMC")& 
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


# Cell chat
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)

table(obj$cell_type)
Idents(obj) <- 'cell_type'
scRNA <- subset(obj, idents = c("sbORS", "UHF", "BG", "TAC", "CP", "bORS", 
                                "SMC", "FB", "EC", "ML",
                                "HMC", "IRS", "CO","FC",
                                "TC","Mac","LC",
                                "BC1","SwG","BC2","GC","MC","SbG","SC"))
scRNA$cell_type <- as.factor(as.character(scRNA$cell_type))
table(scRNA$orig.ident)
Idents(scRNA) <- 'orig.ident'
sco.F <- subset(scRNA, idents = c('F1', 'F2', 'F3'))
sco.O <- subset(scRNA, idents = c('O1', 'O2', 'O3'))
sco.F <- JoinLayers(sco.F)
sco.O <- JoinLayers(sco.O)

a <- sco.F[["RNA"]]@features %>% rownames()
rownames(sco.F@assays$RNA@layers$data) <- a
b <-sco.F[["RNA"]]@cells %>%rownames()
colnames(sco.F@assays$RNA@layers$data)<-b

a <- sco.O[["RNA"]]@features %>% rownames()
rownames(sco.O@assays$RNA@layers$data) <- a
b <-sco.O[["RNA"]]@cells %>%rownames()
colnames(sco.O@assays$RNA@layers$data)<-b

### Create the cellchat objects
cco.F <- createCellChat(sco.F@assays$RNA@layers$data, meta = sco.F@meta.data, group.by = "cell_type")
cco.O <- createCellChat(sco.O@assays$RNA@layers$data, meta = sco.O@meta.data, group.by = "cell_type")
save(cco.til, cco.pbmc, file = "cco.rda")
CellChatDB <- CellChatDB.human
str(CellChatDB)
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$complex)
head(CellChatDB$cofactor)
head(CellChatDB$geneInfo)
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB
cco.F@DB <- CellChatDB.use
cco.F <- subsetData(cco.F)        # subset the expression data of signaling genes for saving computation cost
cco.F <- identifyOverExpressedGenes(cco.F) 
cco.F <- identifyOverExpressedInteractions(cco.F)
cco.F <- projectData(cco.F, PPI.human)                              

#  Assume the possibility of the chating between cells according to the gene expression value
cco.F <- computeCommunProb(cco.F, raw.use = FALSE, population.size = TRUE) 
#  Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cco.F <- filterCommunication(cco.F, min.cells = 10)
F.net <- subsetCommunication(cco.F)
write.csv(F.net, "cell_signaling_net_F.csv")
#Signaling strength assuming
cco.F <- computeCommunProbPathway(cco.F)
F.netlevel <- subsetCommunication(cco.F, slot.name = "netF")
write.csv(F.netlevel, "net_strength_F.csv")

cco.F <- aggregateNet(cco.F)
cco.F <- netAnalysis_computeCentrality(cco.F, slot.name = "netP")
cco.F <- computeNetSimilarity(cco.F, type = "functional")
cco.F <- computeNetSimilarity(cco.F, type = "structural")

cco.O@DB <- CellChatDB.use
cco.O <- subsetData(cco.O)                                                       
cco.O <- identifyOverExpressedGenes(cco.O) 
cco.O <- identifyOverExpressedInteractions(cco.O)
cco.O <- projectData(cco.O, PPI.human)                              
cco.O <- computeCommunProb(cco.O, raw.use = FALSE, population.size = TRUE)
#  Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cco.O <- filterCommunication(cco.O, min.cells = 10)
O.net <- subsetCommunication(cco.O)
write.csv(O.net, "cell_signaling_net_O.csv")
#Signaling strength assuming
cco.O <- computeCommunProbPathway(cco.O)
O.netlevel <- subsetCommunication(cco.O, slot.name = "netF")
write.csv(O.netlevel, "net_strength_O.csv")

cco.O <- aggregateNet(cco.O)
cco.O <- netAnalysis_computeCentrality(cco.O, slot.name = "netP")
cco.O <- computeNetSimilarity(cco.O, type = "functional")
cco.O <- netClustering(cco.O, type = "functional")
cco.O <- computeNetSimilarity(cco.O, type = "structural")

cco.list <- list(Forehead = cco.F, Occipital = cco.O)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat, measure = "weight")
h1+h2

par(mfrow = c(1,1))
weight.max <- getMaxWeight(cco.list, attribute = c("idents","count"))
for (i in 1:length(cco.list)) {
  netVisual_circle(cco.list[[i]]@net$count, weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 6, 
                   title.name = paste0("Number of interactions - ", names(cco.list)[i]))
}

netVisual_circle(cco.list[[2]]@net$count, weight.scale = T, label.edge= F, 
                 edge.weight.max = weight.max[2], edge.width.max = 6, 
                 title.name = paste0("Number of interactions - ", names(cco.list)[2]))


## Comparison of pathways strength
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

install.packages("uwot")
library(uwot)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, umap.method = 'uwot', type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, umap.method = "uwot", type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
rankSimilarity(cellchat, type = "functional") + ggtitle("Functional similarity of pathway")

# Cellular signaling modes comparison
library(ComplexHeatmap)
# Comparison of whole signals
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "all", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 16, height = 30)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "all", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 16, height = 30)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
#  Compare_signal_pattern_outgoing
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "outgoing", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 10, height = 30)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "outgoing", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 10, height = 30)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
#  Compare_signal_pattern_incoming
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 10, height = 30)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "incoming", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 10, height = 30)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))

#Ligand-Receptor comparison analysis
# Cell type level order:BC EC FB GL HFB LC Mac MC ML SbG SC SwG TC UHF
levels(cellchat@idents$joint)
netVisual_bubble(cellchat, sources.use = c(4,5,14), targets.use = "EC",  comparison = c(1, 2), angle.x = 45)
p1 <- netVisual_bubble(cellchat, sources.use =  c(4,5,10,14), targets.use = "EC", comparison = c(1, 2), signaling = c("EPHA","EPHB"),
                       max.dataset = 2, title.name = "Increased signaling to EC", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use =  c(4,5,10,14), targets.use = "EC", comparison = c(1, 2), signaling = c("EPHA","EPHB"),
                       max.dataset = 1, title.name = "Decreased signaling to EC", angle.x = 45, remove.isolate = T)
p1
p2
# Chord plot
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_chord_gene(cco.list[[i]], sources.use = c(4,5,14), targets.use = "EC", slot.name = "netP",
                       lab.cex = 0.1, legend.pos.x = 10, legend.pos.y = 20,
                       title.name = paste0("Signaling to EC - ", names(cco.list)[i]))
}
netVisual_chord_gene(cco.list[[1]], sources.use = c(4,5,14), targets.use = "EC", slot.name = "netP",
                     lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 20,
                     title.name = paste0("Signaling to EC - ", names(cco.list)[1]))

par(mfrow = c(1,2))
s.cell <- c("EC", "UHF", "HFB", "GL")
count1 <- cco.list[[1]]@net$count[s.cell, s.cell]
count2 <- cco.list[[2]]@net$count[s.cell, s.cell]
weight.max <- max(max(count1), max(count2))
netVisual_circle(count1, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Number of interactions-", names(cco.list)[1]))
netVisual_circle(count2, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Number of interactions-", names(cco.list)[2]))

netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")


# Get the subset of Forehead group
cco.F_filtered <- subsetCellChat(cco.F, idents.use = select.cell)
# Get the subset of Occipital group
cco.O_filtered <- subsetCellChat(cco.list$Occipital, idents.use = select.cell)

pathways.show <- c("EPHA") 
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show) 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(cco.list)[i]))
}
ht <- list()
for (i in 1:length(cco.list)) {
  ht[[i]] <- netVisual_heatmap(cco.list[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(cco.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(1.0, "cm"))

num.link <- sapply(cco.list, function(x) {rowSums(x@net$count) +
colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(cco.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(cco.list[[i]],
                                               title = names(cco.list)[i], 
                                               weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "EC", do.label = TRUE, top.label = 1)
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "UHF")
patchwork::wrap_plots(plots = list(gg1,gg2))
gg1
gg2


# Display the chosen pathway
groupSize_Forehead <- as.numeric(table(cellchat@idents$Forehead))
groupSize_Occipital <- as.numeric(table(cellchat@idents$Occipital))
pathways.show <- c("EPHA") 
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[i]))
}
netVisual_aggregate(cco.list[[1]], vertex.weight = groupSize_Forehead, signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[1]))
netVisual_aggregate(cco.list[[2]], vertex.weight = groupSize_Occipital, signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[2]))

# heatmap
pathways.show <- c("VEGF") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(cco.list)) {
  ht[[i]] <- netVisual_heatmap(cco.list[[i]], signaling = pathways.show, 
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(cco.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(1.0, "cm"))

saveRDS(obj,file = "obj.rds")
saveRDS(cellchat,file = "cellchat.rds")
saveRDS(cco.O,file = "cco.O.rds")
saveRDS(cco.F,file = "cco.F.rds")
saveRDS(cco.list,file = "cco.list.rds")

