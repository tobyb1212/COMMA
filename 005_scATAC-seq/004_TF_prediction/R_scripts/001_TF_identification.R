library(ArchR)
library(tidyverse)
library(ComplexHeatmap)
library(chromVAR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)

setwd("/COMMA/results/scatac/transcription_factors")

# Load Project
proj_vascular <- loadArchRProject(path="/COMMA/results/scatac/ArchR/vascular_clustered", showLogo=FALSE)

# Add JASPAR 2026 Motifs
pfm <- TFBSTools::readJASPARMatrix("/COMMA/results/scatac/transcription_factors/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt", matrixClass = "PFM")
pwm <- toPWM(pfm)
names(pwm) <- make.unique(names(pwm), sep="_")
proj_vascular <- addMotifAnnotations(ArchRProj = proj_vascular, motifPWMs = pwm, annoName = "Motif", force = TRUE)

# Define Motifs Enriched in Vascular Cell Types
setwd("/COMMA/results/scatac/transcription_factors/vascular_cell_type")
markersPeaks <- getMarkerFeatures(ArchRProj = proj_vascular, useMatrix = "PeakMatrix", groupBy = "vascular_cell_type",
      bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
motifsUp_vascular <- peakAnnoEnrichment(seMarker = markersPeaks, ArchRProj = proj_vascular, peakAnnotation = "Motif", cutOff = "Pval <= 0.05 & Log2FC >= 0.5")
saveRDS(motifsUp_vascular, "vascular_cell_type_motifsUp.rds")
cmp_names <- colnames(assay(motifsUp_vascular, assayNames(motifsUp_vascular)[1]))
moi <- c()
for (cmp in cmp_names) {motif_stats <- data.frame(motif = rownames(motifsUp_vascular),
    setNames(lapply(assayNames(motifsUp_vascular), function(a) assay(motifsUp_vascular, a)[, cmp]),
      assayNames(motifsUp_vascular)),check.names = FALSE)
  motif_stats <- motif_stats[order(-motif_stats$mlog10Padj), ]
  motif_stats$rank <- seq_len(nrow(motif_stats))
  moi <- c(moi, head(motif_stats$motif, 5))}

# Supplemental Figure 12A
pal <- colorRampPalette(
  c("#BDBDBD", "#FEE8A8", "#FEC44F", "#FE9929", "#D95F0E", "#7F2704")
)(100)
heatmapEM <- try(plotEnrichHeatmap(motifsUp_vascular, n=5, transpose=TRUE, pal = pal))
try({pdf("vascular_cell_type_motifsUp_heatmap.pdf", width = 8, height = 6); ComplexHeatmap::draw(heatmapEM, heatmap_legend_side="bot", annotation_legend_side="bot"); dev.off()}, silent=TRUE)
try({png("vascular_cell_type_motifsUp_heatmap.png", width = 1600, height = 1200, res = 200, type = "cairo"); ComplexHeatmap::draw(heatmapEM, heatmap_legend_side="bot", annotation_legend_side="bot"); dev.off()}, silent=TRUE)

# Motif deviation analysis (proj_vascular)
proj_vascular <- addBgdPeaks(proj_vascular, force=TRUE)
proj_vascular <- addDeviationsMatrix(ArchRProj = proj_vascular, peakAnnotation = "Motif", force = TRUE)
VarDev <- getVarDeviations(proj_vascular, name = "MotifMatrix", plot=FALSE)
write.csv(as.data.frame(VarDev), "MotifDeviation_ranked_vascular.csv", row.names = FALSE)
png("MotifDeviation_ranked_vascular.png", 2000, 1600, res = 200, type = "cairo")
getVarDeviations(proj_vascular, name = "MotifMatrix", plot=TRUE)
dev.off()
markerMotifs_dev <- VarDev$name[1:5]
motifs_vasc_all <- unique(c(moi, markerMotifs_dev))
markerMotifs_vasc <- getFeatures(proj_vascular, select = paste(motifs_vasc_all, collapse = "|"), useMatrix = "MotifMatrix")
saveRDS(markerMotifs_vasc, "markerMotifs_vascular_cell_type.rds")
markerMotifs_vasc_z <- markerMotifs_vasc[grepl("^z:", markerMotifs_vasc)]
outdir <- "deviations"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
proj_vascular$vascular_cell_type_ordered <- dplyr::recode(
    proj_vascular$vascular_cell_type,
    "Arterial Endothelial Cell" = "1_Arterial Endothelial Cell",
    "Capillary Endothelial Cell" = "2_Capillary Endothelial Cell",
    "Venous Endothelial Cell" = "3_Venous Endothelial Cell",
    "Vascular Smooth Muscle Cell" = "Vascular Smooth Muscle Cell",
    "Pericyte" = "5_Pericyte"
)
# Supplemental Figure 12B (SOX17) and Figure 5F (SRF)
for (tf in c("SOX17", "SRF")){
  p <- NULL
  tryCatch({p <- plotGroups(ArchRProj = proj_vascular, groupBy = "vascular_cell_type", colorBy = "MotifMatrix", name = tf, imputeWeights = getImputeWeights(proj_vascular))}, error = function(e) message("Failed: ", tf))
  if (is.null(p)) next
  tryCatch({safe_tf <- gsub("[^A-Za-z0-9_:-]", "_", tf)
    pdf(file.path(outdir, paste0("groups_",safe_tf, ".pdf")), width = 10, height = 8); print(p); dev.off()
    png(file.path(outdir, paste0("groups_",safe_tf, ".png")), width = 2000, height = 1600, res = 200, type = "cairo"); print(p); dev.off()
  }, error = function(e) message("Failed: ", tf))}

# Footprinting (proj_vascular)
motifPositions_v <- getPositions(proj_vascular)
marker_positions_v <- unlist(lapply(motifs_vasc_all, function(x) grep(x, names(motifPositions_v), value = TRUE)))
seFoot_v <- getFootprints(ArchRProj = proj_vascular, positions = motifPositions_v[marker_positions_v], groupBy = "vascular_cell_type")
saveRDS(seFoot_v, "footprinting_vascular_cell_type.rds")
outdir="footprints"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
plots <- NULL
try(plots <- plotFootprints(seFoot = seFoot_v, ArchRProj = proj_vascular, normMethod = "Subtract", plotName = "Footprints_Vascular_Cell_Type_Subtract_Bias", addDOC = FALSE, smoothWindow = 5, plot=FALSE, force=TRUE))
# Supplemental Figure 12D
if (!is.null(plots)) {
  for (tf in "SRF") {
    safe_tf <- gsub("[^A-Za-z0-9_:-]", "_", tf)
    p <- plots[[tf]]
    tryCatch({pdf(file.path(outdir, paste0("footprints_", safe_tf, ".pdf")),width = 10, height = 8); grid.draw(p); dev.off()
      png(file.path(outdir, paste0("footprints_", safe_tf, ".png")), width = 2000, height = 1600, res = 200, type = "cairo"); grid.draw(p); dev.off()}, error = function(e) message("Failed: ", tf))}}