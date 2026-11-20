library(tidyverse)
library(tidyr)
library(rtracklayer)
library(GenomicRanges)
library(ArchR)
library(SummarizedExperiment)
library(data.table)
library(Matrix)

# Use Differentially Accessible Peaks from Vascular Cell Types
files <- list.files("/COMMA/results/scatac/gwas_overlaps/vascular", pattern = "*.csv", full.names = TRUE)
files <- files[files != "/COMMA/results/scatac/gwas_overlaps/vascular/pseudobulk_counts.csv"]
files <- files[files != "/COMMA/results/scatac/gwas_overlaps/vascular/pseudobulk_metadata.csv"]

# Make these values continuous annotations
for (file in files) {
  bed4 <- read.csv(file) %>%
    separate(feature, into = c("chr", "start", "end"), sep = "_", convert = TRUE) %>%
    mutate(p = pmax(Pval, .Machine$double.xmin), score = -log10(p) * Log2FC, chr = paste0("chr", chr)) %>%
    transmute(chr, start = start - 1, end, score)
  new_path <- basename(file) %>%
    str_replace("vascular", "continuous_annotations") %>%
    str_replace("\\.csv$", ".bed")
  new_path <- file.path("/COMMA/results/GWAS/s-LDSC/annotations", new_path)
  write.table(bed4,file = new_path,sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)}