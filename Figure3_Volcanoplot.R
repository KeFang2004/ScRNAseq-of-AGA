# > head(DEG)
#              p_val avg_log2FC pct.1 pct.2     p_val_adj     cluster  gene label
#RPS12 1.273332e-143  0.7298951 1.000 0.991 1.746248e-139 Naive CD4 T RPS12    Up
#RPS6  6.817653e-143  0.6870694 1.000 0.995 9.349729e-139 Naive CD4 T  RPS6    Up
#


color.pals = c("#337ab7", "#4e99c7", "#a7cfe4","#d1e5f0",
               "#21579d", "#2682d1","#549da3","#96cb8f","#4dae47")
#
#' multi volcano plot for scRNA-seq
#' @version 0.2 change legend order
#'
#' @param dat Seurat FindAllMarkers returns, must set only.pos = F;
#' @param color.arr color list, default same as Seurat
#' @param onlyAnnotateUp only annote gene symbols for up genes
#' @param log2Foldchang threshold for annotation
#' @param adjp  threshold for annotation
#' @param top_marker gene number for annotation
#' @param max_overlaps annotation label overlapping
#'
#' @return ggplot2 obj
#' @export
#'
#' @examples
multiVolcanoPlot = function(dat, color.arr=NULL, onlyAnnotateUp=T,
                            log2Foldchang=1.5, adjp=0.05, top_marker=5, 
                            max_overlaps=20, width=0.9){
  library(dplyr)
  library(ggrepel)
  # set default color list
  if(is.null(color.arr)){
    len = length(unique(dat$cluster))
    color.arr=scales::hue_pal()(len)
  }
  
  dat.plot <- dat %>% mutate(
    "significance"=case_when(p_val_adj < adjp & avg_log2FC >= log2Foldchang  ~ 'Up',
                             p_val_adj < adjp & avg_log2FC <= -log2Foldchang  ~ 'Down',
                             TRUE ~ 'None'))
  tbl = table(dat.plot$significance)
  print( tbl )
  background.dat <- data.frame(
    dat.plot %>% group_by(cluster) %>% filter(avg_log2FC>0) %>%
      summarise("y.localup"=max(avg_log2FC)),
    dat.plot %>% group_by(cluster) %>% filter(avg_log2FC<=0) %>%
      summarise("y.localdown"=min(avg_log2FC)),
    x.local=seq(1:length(unique(dat.plot$cluster)))
  ) %>% select(-cluster.1)
  #names(background.dat)
  #head(background.dat)
  #dim(background.dat)
  
  #
  x.number <- background.dat %>% select(cluster, x.local)
  dat.plot <- dat.plot%>% left_join(x.number,by = "cluster")
  #names(dat.plot)
  #head(dat.plot)
  
  #selecting top-up and top-down proteins
  dat.marked.up <- dat.plot %>% filter(significance=="Up") %>%
    group_by(cluster) %>% arrange(-avg_log2FC) %>%
    top_n(top_marker,abs(avg_log2FC))
  dat.marked.down <- dat.plot %>% filter(significance=="Down") %>%
    group_by(cluster) %>% arrange(avg_log2FC) %>%
    top_n(top_marker,abs(avg_log2FC))
  dat.marked <- dat.marked.up %>% bind_rows(dat.marked.down)
  #referring group information data
  dat.infor <- background.dat %>%
    mutate("y.infor"=rep(0,length(cluster)))
  #names(dat.infor)
  #dim(dat.infor)
  #head(dat.infor)
  
  ##plotting:
  #setting color by loading local color schemes
  vol.plot <- ggplot()+
    # background
    geom_col(background.dat,mapping=aes(x.local, y.localup),
             fill="grey80", alpha=0.2, width=0.9, just = 0.5)+
    geom_col(background.dat,mapping=aes(x.local,y.localdown),
             fill="grey80", alpha=0.2, width=0.9, just = 0.5)+
    # point plot
    geom_jitter(dat.plot, mapping=aes(x.local, avg_log2FC, #x= should be number, Not string or factor
                                      color=significance),
                size=0.8, width = 0.4, alpha= 1)+
    scale_color_manual(name="significance", 
                       breaks = c('Up', 'None', 'Down'),
                       values = c("#FF0D58","#F8F8F8", "#3E90E3")) + #set color for: Down None   Up
    geom_tile(dat.infor, mapping=aes(x.local, y.infor), #x axis color box
              height = log2Foldchang*1.3,
              fill = color.arr[1:length(unique(dat.plot$cluster))],
              alpha = 0.5,
              width=width) +
    labs(x=NULL,y="log2 Fold change")+
    geom_text(dat.infor, mapping=aes(x.local,y.infor,label=cluster))+
    # Down is not recommend, not meaningful, hard to explain; so prefer dat.marked.up to dat.marked
    ggrepel::geom_label_repel(data=if(onlyAnnotateUp) dat.marked.up else dat.marked, #gene symbol, of up group default
                              mapping=aes(x=x.local, y=avg_log2FC, label=gene),
                              force = 1, #size=2,
                              #max.overlaps = max_overlaps,
                              label.size = 0, #no border
                              fill="#00000000", #box fill color
                              seed = 233,
                              min.segment.length = 0,
                              force_pull = 1,
                              box.padding = 0.1,
                              segment.linetype = 3,
                              #segment.color = 'black',
                              #segment.alpha = 0.5,
                              #direction = "x", #line direction
                              hjust = 0.5)+
    annotate("text", x=1.5, y=max(background.dat$y.localup)+1,
             label=paste0("|log2FC|>=", log2Foldchang, " & FDR<", adjp))+
    theme_classic(base_size = 12)+
    
    theme(
      axis.title = element_text(size = 13, color = "black"),
      axis.text = element_text(size = 15, color = "black"),
      axis.line.y = element_line(color = "black", size = 0.8),
      #
      axis.line.x = element_blank(), #no x axis line
      axis.ticks.x = element_blank(), #no x axis ticks
      axis.title.x = element_blank(), #
      axis.text.x = element_blank(),
      #
      legend.spacing.x = unit(0.1,'cm'),
      legend.key.width = unit(0.5,'cm'),
      legend.key.height = unit(0.5,'cm'),
      legend.background = element_blank(),
      legend.box = "horizontal",
      legend.position = c(0.13, 0.77),legend.justification = c(1,0)
    )+
    guides( #color = guide_legend( override.aes = list(size=5) ), #legend circle size
      color=guide_legend( override.aes = list(size=5), title="Change")
    )
  #guides(fill=guide_legend(title="Change"))+ #change legend title
  vol.plot
}

DimPlot(OL, label = T)
ol.markers <- FindAllMarkers(OL, min.pct = 0.25, logfc.threshold = 1.5)
ol.markers %>% group_by(cluster) %>% top_n(n = 2, wt = abs(avg_log2FC) )
write.csv(ol.markers,"Outer_layer_cell_markers.csv")
DEG = ol.markers
multiVolcanoPlot(DEG, color.pals,onlyAnnotateUp = F)

