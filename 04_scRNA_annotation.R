DefaultAssay(GARP.int) 
sel.clust = "SCT_snn_res.1"
GARP.int <- SetIdent(GARP.int, value = sel.clust)
table(GARP.int@active.ident) 

#Dim#
DimPlot(GARP.int, reduction = "umap", group.by = "SCT_snn_res.1",label=T,split.by = "orig.ident")
DotPlot(GARP.int,features =c("Foxp3","Lrrc32"),assay="RNA")

#Marker genes#
all.markers <- FindAllMarkers(GARP.int, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
all.markers<- subset(all.markers, all.markers$p_val_adj < 0.05)
all.markers <- dplyr :: arrange(all.markers, desc(all.markers$avg_log2FC))
all.markers <- dplyr :: arrange(all.markers, all.markers$cluster)
library(openxlsx)
write.xlsx(all.markers,"GARP.int_markers_res1.xlsx")

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
p_all_markers <- DotPlot(GARP.int, features = genes_to_check,assay='RNA')  
p_all_markers

#res1#
celltype=data.frame(ClusterID=0:30,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(0,2,9,29,27),2]='IgM+plasma B' 
celltype[celltype$ClusterID %in% c(22),2]='IgG+plasma B' 
celltype[celltype$ClusterID %in% c(18),2]='Th17' 
celltype[celltype$ClusterID %in% c(8,12,24),2]='Tregs'
celltype[celltype$ClusterID %in% c(1),2]='Naive T' 
celltype[celltype$ClusterID %in% c(7),2]='NKT'
celltype[celltype$ClusterID %in% c(11),2]='NK'
celltype[celltype$ClusterID %in% c(14,23),2]='Neutrophils' 
celltype[celltype$ClusterID %in% c(16),2]='Proliferating T'
celltype[celltype$ClusterID %in% c(13),2]='Naive CD4+T' 
celltype[celltype$ClusterID %in% c(10),2]='Naive CD8+T' 
celltype[celltype$ClusterID %in% c(17),2]='Cytotoxic CD8+T'
celltype[celltype$ClusterID %in% c(19),2]='Gd17T'
celltype[celltype$ClusterID %in% c(25),2]='pDCs'
celltype[celltype$ClusterID %in% c(4),2]='Fn1+MoMacs' 
celltype[celltype$ClusterID %in% c(15),2]='C1qc+MoMacs' 
celltype[celltype$ClusterID %in% c(6,5),2]='CD36+MoMacs' 
celltype[celltype$ClusterID %in% c(3,21),2]='ISG+MoMacs' 
celltype[celltype$ClusterID %in% c(20,28),2]='cDCs' 
celltype[celltype$ClusterID %in% c(26),2]='AMs' 
celltype[celltype$ClusterID %in% c(30),2]='Mast cells'


#Annotate#
celltype 
table(celltype$celltype)

GARP.int@meta.data$celltype
GARP.int@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  GARP.int@meta.data[which(GARP.int@meta.data$SCT_snn_res.1 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(GARP.int@meta.data$celltype)

new.cluster.ids <- celltype$celltype
names(new.cluster.ids) <- as.character(celltype$ClusterID) # 按原ClusterID显式命名，避免levels顺序错配。
GARP.int<- RenameIdents(GARP.int, new.cluster.ids)
p<-DimPlot(GARP.int, reduction = "umap", label = T, pt.size = 0.5,label.size = 5)  +theme(legend.position = "none")+ggtitle("")
p


saveRDS(GARP.int, "GARP.int_annotated.rds")
