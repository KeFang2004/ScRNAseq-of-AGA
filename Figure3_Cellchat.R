###   Cell Chat Analyais   ###
# Mapping to old clusters
Idents(obj, cells = colnames(OL)) <- Idents(OL)
obj$celltype <- Idents(obj)
DimPlot(obj,label = TRUE)
obj <- subset(obj, subset = celltype != "Inner layer cluster")

library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)

table(obj$celltype)
Idents(obj) <- 'celltype'
scRNA <- subset(obj, idents = c("bORS-WI-SERPINF1", "bORS-ECM-COL4A1", "bORS-SD-TNFRSF19", "bORS-TD-PPARGC1A",
                                "sbORS-KR-KRT6B", "sbORS-PL-PARD3", "CP", "bTAC",
                                "BG","ML","IMM",
                                "VASC","FIB","MISC"))
scRNA$celltype <- as.factor(as.character(scRNA$celltype))
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

### 创建cellchat对象
cco.F <- createCellChat(sco.F@assays$RNA@layers$data, meta = sco.F@meta.data, group.by = "celltype")
cco.O <- createCellChat(sco.O@assays$RNA@layers$data, meta = sco.O@meta.data, group.by = "celltype")
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
cco.F <- subsetData(cco.F)                                                       # subset the expression data of signaling genes for saving computation cost
cco.F <- identifyOverExpressedGenes(cco.F) 
cco.F <- identifyOverExpressedInteractions(cco.F)
cco.F <- projectData(cco.F, PPI.human)      
#  根据表达值推测细胞互作的概率（cellphonedb是用平均表达值代表互作强度）。
cco.F <- computeCommunProb(cco.F, raw.use = FALSE, population.size = TRUE) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
#  Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cco.F <- filterCommunication(cco.F, min.cells = 10)
F.net <- subsetCommunication(cco.F)
#Signaling strength assuming
cco.F <- computeCommunProbPathway(cco.F)
F.netlevel <- subsetCommunication(cco.F, slot.name = "netF")
write.csv(F.net, "Figure3_cell_signaling_net_F.csv")
write.csv(F.netlevel, "Figure3_net_strength_F.csv")
cco.F <- aggregateNet(cco.F)
cco.F <- netAnalysis_computeCentrality(cco.F, slot.name = "netP")
cco.F <- computeNetSimilarity(cco.F, type = "functional")
cco.F <- computeNetSimilarity(cco.F, type = "structural")

cco.O@DB <- CellChatDB.use
cco.O <- subsetData(cco.O)                                                       
cco.O <- identifyOverExpressedGenes(cco.O) 
cco.O <- identifyOverExpressedInteractions(cco.O)
cco.O <- projectData(cco.O, PPI.human)                             
#  根据表达值推测细胞互作的概率（cellphonedb是用平均表达值代表互作强度）。
cco.O <- computeCommunProb(cco.O, raw.use = FALSE, population.size = TRUE) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
#  Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cco.O <- filterCommunication(cco.O, min.cells = 10)
O.net <- subsetCommunication(cco.O)
#Signaling strength assuming
cco.O <- computeCommunProbPathway(cco.O)
O.netlevel <- subsetCommunication(cco.O, slot.name = "netF")
write.csv(O.net, "figure3_cell_signaling_net_O.csv")
write.csv(O.netlevel, "figure3_net_strength_O.csv")
cco.O <- aggregateNet(cco.O)
cco.O <- netAnalysis_computeCentrality(cco.O, slot.name = "netP")
cco.O <- computeNetSimilarity(cco.O, type = "functional")
cco.O <- computeNetSimilarity(cco.O, type = "structural")


cco.list <- list(Forehead = cco.F, Occipital = cco.O)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)
#barplot
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
#Netplot
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
#Heatmap
par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat, measure = "weight")
h1+h2
#Number of interactions netplot
par(mfrow = c(1,2))
weight.max <- getMaxWeight(cco.list, attribute = c("idents","count"))
for (i in 1:length(cco.list)) {
  netVisual_circle(cco.list[[i]]@net$count, weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 4, 
                   title.name = paste0("Number of interactions - ", names(cco.list)[i]))
}
netVisual_circle(cco.list[[2]]@net$count, weight.scale = T, label.edge= F, 
                 edge.weight.max = weight.max[2], edge.width.max = 6, 
                 title.name = paste0("Number of interactions - ", names(cco.list)[2]))

## Signaling pathway strength comparison
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

library(uwot)
library(futures)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, umap.method = 'uwot', type = "functional")
cellchat <- netClustering(cellchat, type = "functional",do.parallel = FALSE)
#cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, umap.method = "uwot", type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3)
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
                                        title = names(cco.list)[1], width = 10, height = 30,color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "outgoing", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 10, height = 30,color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
#  Compare_signal_pattern_incoming
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 10, height = 30, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "incoming", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 10, height = 30,color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))


#展示两组间各类细胞incoming与outcoming通讯的强度
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

#展示两组间某类细胞输出与输入信号强度的差异
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "IMM", do.label = TRUE, top.label = 1)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "BG", do.label = TRUE, top.label = 1)
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "FIB")
patchwork::wrap_plots(plots = list(gg1,gg2))
gg1
gg2



#Ligand-Receptor comparison analysis
levels(cellchat@idents$joint)
netVisual_bubble(cellchat, sources.use = "FIB", targets.use = c(1,2,3,4,5,6,7,12,13),  comparison = c(1, 2), angle.x = 45)
p1 <- netVisual_bubble(cellchat, sources.use =  c(8,9,10,11,14), targets.use = "BG", comparison = c(1, 2), signaling = c("NOTCH","EPHA"),
                       max.dataset = 2, title.name = "Increased signaling to BG", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use =  c(8,9,10,11,14), targets.use = "BG", comparison = c(1, 2), signaling = c("BMP","COLLAGEN"),
                       max.dataset = 1, title.name = "Decreased signaling to BG", angle.x = 45, remove.isolate = T)
p1
p2

# Chord plot demonstrating signaling pathways betweend selected cell clusters
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_chord_gene(cco.list[[i]], sources.use = "FIB", targets.use = c(1,2,3,4,5,6,7,13,14), slot.name = "netP",
                       lab.cex = 0.8, legend.pos.x = 10, legend.pos.y = 20,thresh = 0.01,
                       title.name = paste0("Signaling to Outer layer - ", names(cco.list)[i]))
}
netVisual_chord_gene(cco.list[[2]], sources.use = c(8,9,10,11,12,15), targets.use = "BG", slot.name = "netP",
                     lab.cex = 0.5, legend.pos.x = 5, legend.pos.y = 10,
                     title.name = paste0("Signaling to BG - ", names(cco.list)[2]))

# Show interaction number between selected cells
par(mfrow = c(1,2))
s.cell <- c("BG", "FIB", "IMM","MISC","ML","VASC")
count1 <- cco.list[[1]]@net$count[s.cell, s.cell]
count2 <- cco.list[[2]]@net$count[s.cell, s.cell]
weight.max <- max(max(count1), max(count2))
netVisual_circle(count1, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Number of interactions-", names(cco.list)[1]))
netVisual_circle(count2, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Number of interactions-", names(cco.list)[2]))

# Show selected signaling pathway by network
par(mfrow = c(1,2), xpd=TRUE)
groupSize_Forehead <- as.numeric(table(cellchat@idents$Forehead))
groupSize_Occipital <- as.numeric(table(cellchat@idents$Occipital))
pathways.show <- c("AGRN") 
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], 
                      edge.width.max = 10, signaling.name = paste(pathways.show, names(cco.list)[i]))
}


#heatmap
pathways.show <- c("CD45") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(cco.list)) {
  ht[[i]] <- netVisual_heatmap(cco.list[[i]], signaling = pathways.show, 
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(cco.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

### 基因表达情况
cellchat@meta$datasets = factor(cellchat@meta$datasets, 
                                levels = c("Forehead", "Occipital")#设定组别的level
)
plotGeneExpression(cellchat, signaling = "ncWNT", split.by = "datasets", colors.ggplot = T)

