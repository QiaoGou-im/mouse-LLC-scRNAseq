
GARP_merge.filt<-readRDS("GARP_merge.filt_qc.rds")
#split#
table(GARP_merge.filt@meta.data$orig.ident)

GARP.list <- SplitObject(GARP_merge.filt, split.by = "orig.ident")
GARP.list

#SCTranform#
GARP.list <- lapply(X = GARP.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = GARP.list, nfeatures = 3000)
GARP.list <- PrepSCTIntegration(object.list = GARP.list, anchor.features = features)

immune.anchors <- FindIntegrationAnchors(object.list = GARP.list, normalization.method = "SCT",anchor.features = features)
GARP.int <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

GARP.int = RunPCA(GARP.int, npcs = 30)
GARP.int = RunTSNE(GARP.int, npcs = 30)
GARP.int = RunUMAP(GARP.int, dims = 1:30)

saveRDS(GARP.int,"GARP.int_scale.rds")

#Clustering#
GARP.int=FindNeighbors(GARP.int, dims = 1:20, k.param = 60, prune.SNN = 1/15)

for (res in c(0.5,1)) {
  GARP.int=FindClusters(GARP.int, graph.name = "integrated_snn", resolution = res, algorithm = 1)
}
apply(GARP.int@meta.data[,grep("integrated_snn_res",colnames(GARP.int@meta.data))],2,table)

p1_dim=plot_grid(ncol = 3, DimPlot(GARP.int, reduction = "umap", group.by = "integrated_snn_res.0.5",label = T) + 
                   ggtitle("louvain_0.5"), DimPlot(GARP.int, reduction = "umap", group.by = "integrated_snn_res.1",label = T) + 
                   ggtitle("louvain_1"))

p1_dim

ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_low.png",width = 14)

saveRDS(GARP.int, "GARP.int_SCT_annotated.rds")
