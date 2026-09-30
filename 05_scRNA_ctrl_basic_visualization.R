
library(Seurat)
library(ggplot2)
library(patchwork)
GARP.int<-readRDS("GARP.int_annotated_res1.rds")

palette <- c("#D2B5E2", "#84AB53", "#DDC0BD", "#5C5FD7", "#DD8B9E", "#E8E3AD", "#8AEA55", "#6AE3D8", "#6CA2E1", "#DBEB87", "#748B74",
                      "#7F37E8", "#D4995E", "#7EBFD6", "#B05BD8", "#DD544C", "#9887DB", "#E691E5", "#A9E5B7", "#5EE999", "#A34A84", "#CDE5E3",
                      "#E4DD42", "#DE3FDD", "#7B7C9B", "#E64EAC")

cell<-c("CD36+MonoMacs", "Fn1+MonoMacs", "NK", "NKT", "Neutrophils", "Tregs", "ISG+MonoMacs", "IgM+plasma B", "Cytotoxic CD8+T", "Gd17T", "cDCs", "pDCs", "Th17", "Proliferating T", "Naive T", "Naive CD8+T", "Naive CD4+T", "Mast cells", "IgG+plasma B", "C1qc+MonoMacs", "AMs")
named_palette <- setNames(palette[seq_along(cell)], cell)

GARP.list <- SplitObject(GARP.int, split.by = "orig.ident")
Ctrl<-GARP.list$Ctrl
GARP_cko<-GARP.list$`GARP_cKO`


Idents(Ctrl)<-"celltype"
p1<-DimPlot(Ctrl, label = T, pt.size = 0.5, repel = T,cols = named_palette,label.size = 8)+NoLegend()+labs(x = "UMAP1", y = "UMAP2") +theme(text = element_text(size = 32))+ theme(axis.text.x=element_text(size=28))+theme(axis.text.y=element_text(size=28))
p1
ggsave(plot=p1, filename="Ctrl_UMAP.pdf",width = 10,height = 10)


genes_to_check = c('Ptprc', 'Cd3d', 'Cd3e', 'Cd4','Cd8a','Foxp3','Lrrc32','Ifng','Il4','Il17a','Bcl6',#T cells#
                   'Cd19', 'Cd79a', 'Ms4a1',#B cells#
                   'Cd68', 'Cd163', 'Cd14', #Monocytes#
                   'Adgre1','Siglecf','Cx3cr1', #macrophages#
                   'Itgax','H2-Ab1',#dendritic cells#
                   'Kit','Fcer1a',  ### mast cells
                   'Ly6g' , 'Itgam', #Neutrophils,
                   'Ncr1','Klrb1', #NK cells#,
                   'Mki67','Apoe','Sell','Cd44','Il7r',"Trbc2"
)
p_all_markers <- DotPlot(Ctrl, features = genes_to_check,assay='RNA')  
p_all_markers

ggsave(plot=p_all_markers, filename="Ctrl_annotation_markers.pdf",width = 16,height = 10)

