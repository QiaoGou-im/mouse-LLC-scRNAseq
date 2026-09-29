rm(list = ls())


library(Seurat)
library(patchwork)
library(ggplot2)
library(SeuratData)
library(dplyr)
library(clustree)
library(cowplot)
library(stringr)  


WorkingPath <- "PATH_TO_YOUR_ANALYSIS_DIRECTORY"
setwd(WorkingPath)
getwd()

GARPctl_1<-Read10X(data.dir = "data/ctl_1")
GARPctl_2<-Read10X(data.dir = "data/ctl_2")

GARPcko_1<-Read10X(data.dir = "data/G_cko_1")
GARPcko_2<-Read10X(data.dir = "data/G_cko_2")

GARPctl_1<- CreateSeuratObject(counts = GARPctl_1, project = "Ctrl", min.cells = 3, min.features = 200)
GARPctl_2<- CreateSeuratObject(counts = GARPctl_2, project = "Ctrl", min.cells = 3, min.features = 200)

GARPcko_1<- CreateSeuratObject(counts = GARPcko_1, project = "GARP_cKO", min.cells = 3, min.features = 200)
GARPcko_2<- CreateSeuratObject(counts = GARPcko_2, project = "GARP_cKO", min.cells = 3, min.features = 200)

GARP_merge = merge(GARPctl_1, y = c(GARPctl_2,GARPcko_1,GARPcko_2), add.cell.ids = c("Ctrl_1","Ctrl_2", "GARP_cKO_1","GARP_cKO_2"),
                   project = "LLC", merge.data = TRUE)

as.data.frame(GARP_merge@assays$RNA@counts[1:10, 1:2])
table(GARP_merge@meta.data$orig.ident) 
head(GARP_merge@meta.data)

#QC#
mito_genes=rownames(GARP_merge)[grep("^mt-", rownames(GARP_merge))] 
mito_genes
GARP_merge=PercentageFeatureSet(GARP_merge, "^mt-", col.name = "percent_mito")
fivenum(GARP_merge@meta.data$percent_mito)

ribo_genes=rownames(GARP_merge)[grep("^Rp[sl]", rownames(GARP_merge),ignore.case = T)]
ribo_genes
GARP_merge=PercentageFeatureSet(GARP_merge, "^Rp[sl]", col.name = "percent_ribo")
fivenum(GARP_merge@meta.data$percent_ribo)

rownames(GARP_merge)[grep("^Hb[^(p)]", rownames(GARP_merge),ignore.case = T)]
GARP_merge=PercentageFeatureSet(GARP_merge, "^Hb[^(p)]", col.name = "percent_hb")
fivenum(GARP_merge@meta.data$percent_hb)

feats <- c("nFeature_RNA", "nCount_RNA","percent_mito")
p1<-VlnPlot(GARP_merge, group.by = "orig.ident", features = feats, pt.size = 0.01, ncol = 2) + 
  NoLegend()
p1

feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2=VlnPlot(GARP_merge, group.by = "orig.ident", features = feats, pt.size = 0.01, ncol = 3, same.y.lims=T) + 
  scale_y_continuous(breaks=seq(0, 100, 5)) +
  NoLegend()
p2 

p3=FeatureScatter(GARP_merge, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
p3

selected_c <- WhichCells(GARP_merge, expression = nFeature_RNA > 300)
selected_f <- rownames(GARP_merge)[Matrix::rowSums(GARP_merge@assays$RNA@counts > 0 ) > 3]
GARP_merge.filt <- subset(GARP_merge, features = selected_f, cells = selected_c)
dim(GARP_merge) #19161 34029#
dim(GARP_merge.filt) #18511 33838#

selected_mito <- WhichCells(GARP_merge.filt, expression = percent_mito < 20)
selected_ribo <- WhichCells(GARP_merge.filt, expression = percent_ribo < 60)
selected_hb <- WhichCells(GARP_merge.filt, expression = percent_hb < 1)
length(selected_hb)  #33834#
length(selected_ribo) #33838#
length(selected_mito) #32664#

GARP_merge.filt<- subset(GARP_merge.filt, cells = selected_mito)

GARP_merge.filt <- subset(GARP_merge.filt, cells = intersect(colnames(GARP_merge.filt), selected_ribo))
dim(GARP_merge.filt) # 18511 33125#
table(GARP_merge.filt$orig.ident) #Ctrl  17360 , GARP cKO  15304  # 

#Filtrate the MALAT1 gene,mt- genes, and !Hsp genes#
GARP_merge.filt <- GARP_merge.filt[!grepl("Malat1", rownames(GARP_merge.filt),ignore.case = T), ]
GARP_merge.filt <- GARP_merge.filt[!grepl("^mt-", rownames(GARP_merge.filt),ignore.case = T), ]#
GARP_merge.filt <- GARP_merge.filt[!grepl("Hsp", rownames(GARP_merge.filt),ignore.case = T), ]
GARP_merge.filt <- GARP_merge.filt[!grepl("^Rp([0-9]+-|[ls])", rownames(GARP_merge.filt),ignore.case = T), ]#
dim(GARP_merge.filt) #dim  18367 32664#

#Find viraible features#
GARP_merge.filt = FindVariableFeatures(GARP_merge.filt)
top10 <- head(VariableFeatures(GARP_merge.filt), 10)
plot1 <- VariableFeaturePlot(GARP_merge.filt)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
