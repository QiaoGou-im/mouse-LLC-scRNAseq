
library(Seurat)
library(ggplot2)
library(patchwork)
GARP.int<-readRDS("GARP.int_annotated_res1.rds")

palette <- c("#D2B5E2", "#84AB53", "#DDC0BD", "#5C5FD7", "#DD8B9E", "#E8E3AD", "#8AEA55", "#6AE3D8", "#6CA2E1", "#DBEB87", "#748B74",
                      "#7F37E8", "#D4995E", "#7EBFD6", "#B05BD8", "#DD544C", "#9887DB", "#E691E5", "#A9E5B7", "#5EE999", "#A34A84", "#CDE5E3",
                      "#E4DD42", "#DE3FDD", "#7B7C9B", "#E64EAC")

cell<-c("CD36+MoMacs", "Fn1+MoMacs", "NK", "NKT", "Neutrophils", "Tregs", "ISG+MoMacs", "IgM+plasma B", "Cytotoxic CD8+T", "Gd17T", "cDCs", "pDCs", "Th17", "Proliferating T", "Naive T", "Naive CD8+T", "Naive CD4+T", "Mast cells", "IgG+plasma B", "C1qc+MoMacs", "AMs")
named_palette <- setNames(palette[seq_along(cell)], cell)

DimPlot(GARP.int, label = T, pt.size = 0.5, repel = T,cols = named_palette,split.by = "orig.ident",label.size = 8)+NoLegend()+labs(x = "UMAP1", y = "UMAP2") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank(),text = element_text(size=24))+theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) 
ggsave("Umap2_res1.pdf",width = 20,height = 10)

Idents(GARP.int)<-"celltype"
p1<-DimPlot(GARP.int, label = T, pt.size = 0.5, repel = T,cols = named_palette,label.size = 8)+NoLegend()+labs(x = "UMAP1", y = "UMAP2") +theme(text = element_text(size = 32))+ theme(axis.text.x=element_text(size=28))+theme(axis.text.y=element_text(size=28))
p1
Idents(GARP.int)<-"orig.ident"
p2<-DimPlot(GARP.int, pt.size = 0.5, repel = T,label.size = 8,cols =c("#4DBBD5B3", "#E64B35B3"))+labs(x = "UMAP1", y = "UMAP2") +theme(text = element_text(size = 32))+theme(axis.text.x=element_text(size=28))+theme(axis.text.y=element_text(size=28))
p2
p1+p2      
ggsave(filename = "UMAP_JSI.pdf",width = 10, height = 10)
ggsave(filename = "UMAP_LLC_2.pdf",width = 20, height = 10)

Idents(GARP.int)<-"celltype"

library(dplyr)
library(ggstatsplot)
library(gginnards)
PropPlot <- function(object, groupBy, celltype_order = NULL) {
  # (1)获取绘图数据
  plot_data = object@meta.data %>% 
    dplyr::select(orig.ident, dplyr::all_of(groupBy)) %>%
    dplyr::rename(group = dplyr::all_of(groupBy)) # 修正字符串列名选择；不改数据内容
  
  # Reorder celltype if custom order is provided
  if (!is.null(celltype_order)) {
    plot_data$group <- factor(plot_data$group, levels = celltype_order)
  }
  
  # (2)绘图
  figure = ggbarstats(data = plot_data, 
                      x = group, y = orig.ident,
                      results.subtitle = FALSE,
                      bf.message = FALSE,
                      proportion.test = FALSE,
                      label.args = list(size = 3, 
                                        fill = 'white', 
                                        alpha = 0.85,
                                        fontface = 'bold'),
                      perc.k = 2,
                      title = '',
                      xlab = '',
                      legend.title = 'Celltype',
                      ggtheme = ggpubr::theme_pubclean()) +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(color = 'black', lineend = 'round'),
          legend.position = 'right',
          axis.text.x = element_text(size = 24, color = 'black'),
          axis.text.y = element_text(size = 24, color = 'black'),
          legend.text = element_text(size = 24, color = 'black'),
          legend.title = element_text(size = 24, color = 'black')) +scale_fill_manual(values=palette)
  
  # (3)去除柱子下面的样本量标识：
  gginnards::delete_layers(x = figure, match_type = 'GeomText')
}

p<-PropPlot(GARP.int, groupBy = 'celltype', celltype_order = c("AMs", "C1qc+MoMacs","IgG+plasma B", "Mast cells", "Naive CD4+T", "Naive CD8+T", "Naive T", "Proliferating T", "Th17", "pDCs", "cDCs", "Gd17T", "Cytotoxic CD8+T", "IgM+plasma B",   "ISG+MoMacs","Tregs", "Neutrophils","NKT","NK", "Fn1+MoMacs","CD36+MoMacs"))

p1<-delete_layers(p, "GeomLabel")
p1

ggsave(plot=p1, filename="Ctrl_vs_cKO_cell_composition.pdf",width = 14,height = 10)

#DEG analysis#
GARP.int$celltype.orig.ident <- paste(Idents(GARP.int), GARP.int$orig.ident, sep = "_")
GARP.int$celltype <- Idents(GARP.int)
Idents(GARP.int) <- "celltype.orig.ident"

DefaultAssay(GARP.int)  <- "RNA"
#options(future.globals.maxSize = 8000 * 1024^2)
GARP.int_RNA <- NormalizeData(GARP.int, verbose = FALSE)

CD36MoMacs.DEG <- FindMarkers(GARP.int_RNA, ident.1 = "CD36+MoMacs_GARP_cKO", ident.2 = "CD36+MoMacs_Ctrl",test.use = "MAST",logfc.threshold = 0,
                        min.pct = 0)
head(CD36MoMacs.DEG, n = 15)
CD36MoMacs.DEG$gene<-rownames(CD36MoMacs.DEG)
library(openxlsx)
write.xlsx(CD36MoMacs.DEG,file="CD36_DEG.xlsx")

