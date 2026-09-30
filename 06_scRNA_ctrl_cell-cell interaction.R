
# Create the CellChat object from normalized RNA expression and the annotated
# Seurat identities.
data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data")
labels <- Idents(seurat_object)
meta <- data.frame(group = labels, row.names = names(labels))

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")

levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))
groupSize

# The original script ultimately used the complete mouse CellChat database.
CellChatDB <- CellChatDB.mouse
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)

# Parameters retained exactly from the original analysis.
cellchat <- computeCommunProb(
  cellchat,
  type = "truncatedMean",
  trim = 0.1,
  population.size = TRUE,
  raw.use = FALSE
)
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat, signaling = c("TGFb"))
write.csv(df.net, "CellChat_TGFb_interactions_Ctrl.csv", row.names = FALSE)

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

pdf("CellChat_global_network_Ctrl.pdf", width = 14, height = 7)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of interactions"
)
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction weights/strength"
)
dev.off()

mat <- cellchat@net$weight
pdf("CellChat_outgoing_networks_by_cell_group_Ctrl.pdf", width = 18, height = 14)
par(mfrow = c(4, 5), xpd = TRUE)
for (i in 7:nrow(mat)) {
  mat2 <- matrix(
    0,
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )
  mat2[i, ] <- mat[i, ]
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    edge.weight.max = max(mat),
    title.name = rownames(mat)[i]
  )
}
dev.off()

pathways.show <- c("TGFb")

pdf("Cellchat_Tgfb_Ctrl.pdf", width = 10, height = 10)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()

pdf("CellChat_TGFb_contribution_Ctrl.pdf", width = 8, height = 6)
print(netAnalysis_contribution(cellchat, signaling = pathways.show))
dev.off()

pairLR.Tgf <- extractEnrichedLR(
  cellchat,
  signaling = pathways.show,
  geneLR.return = FALSE
)
LR.show <- pairLR.Tgf[1, ]

pdf("CellChat_TGFb_top_LR_pair_Ctrl.pdf", width = 10, height = 10)
netVisual_individual(
  cellchat,
  signaling = pathways.show,
  pairLR.use = LR.show,
  layout = "chord"
)
dev.off()

# Numeric source and target indices are retained from the original script.
pdf("CellChat_TGFb_selected_cell_network_Ctrl.pdf", width = 12, height = 10)
netVisual_chord_cell(
  cellchat,
  sources.use = 13,
  targets.use = c(1:24),
  lab.cex = 1,
  legend.pos.y = 30,
  signaling = pathways.show
)
dev.off()

pdf("CellChat_selected_gene_network_Ctrl.pdf", width = 10, height = 8)
netVisual_chord_gene(
  cellchat,
  sources.use = 7,
  targets.use = c(5),
  lab.cex = 2,
  legend.pos.y = 30,
  reduce = "0.015"
)
dev.off()

pdf("CellChat_selected_bubble_plots_Ctrl.pdf", width = 12, height = 8)
p1 <- netVisual_bubble(
  cellchat,
  sources.use = c(1:25),
  targets.use = c(21),
  remove.isolate = FALSE,
  direction = -1,
  font.size = 18,
  font.size.title = 18,
  max.dataset = 2
)
print(p1)

p2 <- netVisual_bubble(
  cellchat,
  sources.use = c(7),
  targets.use = c(5),
  remove.isolate = TRUE
)
print(p2)

p3 <- netVisual_bubble(
  cellchat,
  sources.use = c(5),
  targets.use = c(4),
  remove.isolate = FALSE
)
print(p3)
dev.off()

pdf("CellChat_TGFb_receptor_expression_Ctrl.pdf", width = 10, height = 7)
print(
  plotGeneExpression(
    cellchat,
    signaling = "TGFb",
    features = c("Tgfbr1", "Tgfbr2")
  )
)
dev.off()

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

pdf("CellChat_TGFb_signaling_roles_Ctrl.pdf", width = 12, height = 8)
netAnalysis_signalingRole_network(
  cellchat,
  signaling = pathways.show,
  width = 12,
  height = 2.5,
  font.size = 10
)
dev.off()

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht <- netAnalysis_signalingRole_heatmap(
  cellchat,
  signaling = c("TGFb"),
  pattern = "incoming"
)

pdf("CellChat_signaling_role_heatmaps_Ctrl.pdf", width = 14, height = 8)
ComplexHeatmap::draw(ht1 + ht2)
ComplexHeatmap::draw(ht)
dev.off()

gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("TGFb"))

pdf("CellChat_signaling_role_scatter_Ctrl.pdf", width = 14, height = 7)
print(gg1 + gg2)
dev.off()

pdf("CellChat_outgoing_pattern_number_Ctrl.pdf", width = 8, height = 6)
selectK(cellchat, pattern = "outgoing")
dev.off()

nPatterns <- 3
cellchat <- identifyCommunicationPatterns(
  cellchat,
  pattern = "outgoing",
  k = nPatterns
)

pdf("CellChat_outgoing_patterns_Ctrl.pdf", width = 12, height = 8)
print(netAnalysis_dot(cellchat, pattern = "outgoing"))
print(netAnalysis_river(cellchat, pattern = "outgoing"))
dev.off()

pdf("CellChat_TGFb_sender_receiver_Ctrl.pdf", width = 10, height = 8)
netAnalysis_signalingRole_network(
  cellchat,
  signaling = pathways.show,
  width = 8,
  height = 2.5,
  font.size = 10,
  measure = c("outdeg", "indeg"),
  measure.name = c("Sender", "Receiver")
)
dev.off()

cell.chat.ctrl <- cellchat
saveRDS(cell.chat.ctrl, "cellchat_ctrl.rds")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

message("Control CellChat analysis completed. Output directory: ", getwd())
