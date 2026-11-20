library(data.table)
library(tidyverse)
library(qqman)
library(ggrastr)
library(ggrepel)
library(ggplot2)

setwd("/COMMA/results/GWAS/CMVD/finemapping")

manhattan_plot <- function(df, p_col = "p", filename = "manhattan") {
  d <- as.data.table(df)
  # Extract SNP Positions
  tmp <- tstrsplit(d$SNP, ":", keep = 1:2)
  d[, CHR := sub("(?i)^chr", "", tmp[[1]], perl = TRUE)]
  d[, BP := as.integer(tmp[[2]])]
  # Order Chromosomes
  chr_order <- c(as.character(1:22), "X", "Y", "MT")
  d[, CHR := factor(CHR, levels = chr_order, ordered = TRUE)]
  # Compute Cumulative Positions
  setorder(d, CHR, BP)
  chr_sizes <- d[, .(maxBP = max(BP, na.rm = TRUE)), by = CHR]
  chr_sizes[, offset := shift(cumsum(as.numeric(maxBP)), fill = 0)]
  d <- merge(d, chr_sizes[, .(CHR, offset)], by = "CHR", sort = FALSE)
  d[, BPcum := BP + offset]
  # Compute P Values
  if (!p_col %in% names(d)) stop("p_col not found in dataframe")
  d[, Pval := as.numeric(get(p_col))]
  d <- d[!is.na(Pval) & !is.na(BPcum)]
  d[, neglogp := -log10(Pval)]
  d[, is_lead := as.logical(is_lead)]
  axis_set <- d[, .(center = (min(BPcum) + max(BPcum)) / 2), by = CHR]
  # Build Manhattan Plot
  p <- ggplot(d, aes(x = BPcum, y = neglogp, colour = factor(as.integer(as.factor(CHR)) %% 2))) +
    ggrastr::rasterise(geom_point(size = 0.6, alpha = 0.8), dpi = 600) +
    ggrastr::rasterise(geom_point(data = d[is_lead == TRUE], shape = 18, colour = "deeppink2", size = 3), dpi = 600) +
    ggrepel::geom_text_repel(data = d[is_lead == TRUE], aes(label = paste0(CHR, ":", BP)),
                             colour = "deeppink2", size = 2.5, max.overlaps = Inf, show.legend = FALSE) +
    scale_colour_manual(values = c("black", "grey60"), guide = "none") +
    scale_x_continuous(labels = axis_set$CHR, breaks = axis_set$center) + theme_classic() +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", colour = "red") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(x = "Chromosome", y = "-log10(p)") +
    scale_y_continuous(breaks = seq(0, ceiling(max(d$neglogp)), by = 2), limits = c(0, NA))
  # Save File
  outpath <- paste0(filename, ".pdf")
  ggsave(filename = outpath, plot = p, width = 8, height = 6, device = cairo_pdf)
  message("Saved: ", outpath)
}

leads <- read.csv("susie_lead_variants.tsv", sep="\t")
sumstats <- read.table("/COMMA/results/GWAS/CMVD/sumstats/NTCV_Supplementary_Table4.csv",header=TRUE, sep="", stringsAsFactors=FALSE)
# Mark lead variants
leads$SNP <- sub("^chr", "", leads$SNP)
sumstats$is_lead <- sub("^chr", "", sub(":[^:]+:[^:]+$", "", sumstats$SNP)) %in% leads$SNP

# Supplementary Figure 14B
manhattan_plot(sumstats, "P", "sumstats_manhattan")