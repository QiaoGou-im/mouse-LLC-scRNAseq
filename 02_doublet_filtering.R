
#Doublet filtration#
options(future.globals.maxSize = 8000 * 1024^2)
GARP_merge.filt = NormalizeData(GARP_merge.filt)
GARP_merge.filt = ScaleData(GARP_merge.filt, 
                            vars.to.regress = c("nFeature_RNA", "percent_mito"))
GARP_merge.filt = RunPCA(GARP_merge.filt, npcs = 30)
GARP_merge.filt = RunTSNE(GARP_merge.filt, npcs = 30)
GARP_merge.filt = RunUMAP(GARP_merge.filt, dims = 1:30)


nExp <- round(ncol(GARP_merge.filt) * 0.04) 
library(DoubletFinder)
GARP_merge.filt <- doubletFinder_v3(GARP_merge.filt, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
DF.name = colnames(GARP_merge.filt@meta.data)[grepl("DF.classification", colnames(GARP_merge.filt@meta.data))]
p5.dimplot=cowplot::plot_grid(ncol = 2, DimPlot(GARP_merge.filt, group.by = "orig.ident") + NoAxes(), 
                              DimPlot(GARP_merge.filt, group.by = DF.name) + NoAxes())
p5.dimplot
ggsave(filename="doublet_dimplot.pdf",plot=p5.dimplot)
p5.vlnplot=VlnPlot(GARP_merge.filt, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
p5.vlnplot
ggsave(filename="doublet_vlnplot.pdf",plot=p5.vlnplot)

#Filter doublet#
GARP_merge.filt=GARP_merge.filt[, GARP_merge.filt@meta.data[, DF.name] == "Singlet"]
dim(GARP_merge.filt) #18367 31357#
dim(GARP_merge) #19161 34029#

saveRDS(GARP_merge.filt, "GARP_merge.filt_qc.rds")
