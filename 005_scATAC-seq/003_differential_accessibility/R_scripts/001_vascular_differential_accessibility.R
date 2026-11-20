library(ArchR)
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(SummarizedExperiment)
library(GenomicRanges)
library(edgeR)
library(data.table)
library(Matrix)

run_archr_wilcoxon_peaks <- function(proj_path, group_by, outdir, useMatrix = "PeakMatrix", bias = c("TSSEnrichment", "log10(nFrags)"),
  cutoff_pval = 0.05, cutoff_logfc = 0.5, returnGR = TRUE, verbose = TRUE) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  proj <- loadArchRProject(path = proj_path)
  proj@cellColData$Cardiac_Disease <- proj@cellColData$`Cardiac Disease`
  seMarker <- getMarkerFeatures(ArchRProj = proj, groupBy = group_by, 
    useMatrix = useMatrix, bias = bias, testMethod = "wilcoxon", verbose = FALSE)
  saveRDS(seMarker, file.path(outdir, paste0(group_by, "_Marker.rds")))
  cut_off <- paste0("Pval <= ", cutoff_pval, " & Log2FC >= ", cutoff_logfc)
  markers <- getMarkers(seMarker, cutOff = cut_off, returnGR = returnGR)
  saveRDS(markers, file.path(outdir, paste0(group_by, "_markers.rds")))
  # Write one BED per group if PeakMatrix + returnGR
  if (identical(useMatrix, "PeakMatrix") && isTRUE(returnGR)) {
    all_peaks <- getMarkers(seMarker, cutOff = "Pval <= 1",returnGR = TRUE)[[1]]
    sig_peaks <- c()
    for (nm in names(markers)) {
      grp_gr <- unique(markers[[nm]])
      sig_peaks <- c(sig_peaks, grp_gr)
      export(grp_gr, file.path(outdir, paste0(nm, ".bed")), format = "BED")
      grp_all <- getMarkers(seMarker, cutOff = "Pval <= 1",returnGR = TRUE)[nm]
      rd <- as.data.frame(rowData(seMarker))
      pval <- assays(seMarker)[["Pval"]][, nm]
      lfc  <- assays(seMarker)[["Log2FC"]][, nm]
      mdif <- assays(seMarker)[["MeanDiff"]][, nm]
      keep <- pval <= 0.05 & lfc > 0.5
      df <- data.frame(
        feature = paste0(gsub("chr", "",rd$seqnames[keep]), "_", rd$start[keep], "_", rd$end[keep]),
        Pval      = pval[keep],
        Log2FC    = lfc[keep],
        MeanDiff  = mdif[keep],
        stringsAsFactors = FALSE)
      write.csv(df, file.path(outdir, paste0(nm, ".csv")), row.names = FALSE)
    }
    sig_union <- GenomicRanges::reduce(do.call(c, sig_peaks))
    bg_peaks <- GenomicRanges::setdiff(all_peaks, sig_union, ignore.strand = TRUE)
    export(bg_peaks, file.path(outdir, "background.bed"), format = "BED")
  }
}

# Identify Differential Peaks for Vascular Cell Types
# Supplemental Table 10
run_archr_wilcoxon_peaks(proj_path = '/COMMA/results/scatac/ArchR/vascular_clustered',group_by = 'vascular_cell_type', 
  outdir = '/COMMA/results/scatac/gwas_overlaps/vascular/')
