#!/usr/bin/env Rscript
library(tidyverse)
library(edgeR)
library(limma)
library(RRHO)
library(lattice)
library(data.table)

# Define Directory & Conditions to Compare
base_dir <- "/COMMA/results/scrna/downstream_analysis/DGE/pseudobulk/" 
counts_fp <- "/COMMA/results/scrna/downstream_analysis/DGE/pseudobulk/vsmc/vsmc__Primary_Cardiac_Disease_byPatient_counts.csv.gz"
condition <- "Primary Cardiac Disease"

# Helper Function
safe_name <- function(x) gsub("[^A-Za-z0-9_.]", "_", x)

# Find Matching Metadata
dirn <- dirname(counts_fp)
fname_base <- sub("_counts\\.csv(?:\\.gz)?$", "", basename(counts_fp))
metadata_fp <- file.path(dirn, paste0(fname_base, "_metadata.csv"))
if (!file.exists(metadata_fp)) {
  message(" - Metadata not found for ", counts_fp, " (expected ", metadata_fp, ") — skipping.")}
message(" - Found metadata: ", metadata_fp)

# Read In Counts
if (grepl("\\.gz$", counts_fp)) {
    counts_dt <- fread(cmd = paste("zcat", shQuote(counts_fp)), data.table = FALSE)
  } else {
    counts_dt <- fread(counts_fp, data.table = FALSE)
  }
  gene_col <- colnames(counts_dt)[1]
  rownames(counts_dt) <- counts_dt[[gene_col]]
  counts_dt[[gene_col]] <- NULL
  counts_mat <- as.matrix(counts_dt)
  mode(counts_mat) <- "numeric"

  # Align Metadata to Counts
  meta <- read.csv(metadata_fp, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(colnames(counts_mat) %in% rownames(meta))) {
    stop("Sample names in counts and metadata don't match for file: ", counts_fp)
  }
  meta <- meta[colnames(counts_mat), , drop = FALSE]

  # Detect Condition Being Compared
  possible <- intersect(colnames(meta), condition)
  if (length(possible) > 0) {
    candidate_in_fname <- possible[sapply(possible,
      function(x) grepl(gsub(" ", "_", x), fname_base, ignore.case = TRUE))]
    if (length(candidate_in_fname) >= 1) {
      group_col <- candidate_in_fname[1]
    } else {
      group_col <- possible[1]
    }
    message(" - Using detected condition column: ", group_col)
}

# Prepare Group Variable
grp_raw <- factor(as.character(meta$`Primary Cardiac Disease`))
levels_grp <- levels(as.factor(meta$`Primary Cardiac Disease`))
level_counts <- table(meta$`Primary Cardiac Disease`)

# Iterate over Groups
if (!"Control" %in% levels_grp) {
  message("Control not found in meta$`Primary Cardiac Disease`; skipping one-vs-Control tests.")
  stop("Primary Cardiac Disease has no Control donors; cannot calculate disease similarity.")
}

control_colname <- paste0("grp", make.names("Control"))
target_levels <- setdiff(levels_grp, "Control")

dge_lists <- list()
for (lvl in target_levels) {
  message("Running Control vs ", lvl)
  # Keep only Control + current level
  keep_samples <- grp_raw %in% c("Control", lvl)
  dat_meta <- meta[keep_samples, , drop = FALSE]
  dat_counts <- counts_mat[, keep_samples, drop = FALSE]
  dat_meta$grp <- factor(as.character(dat_meta[[group_col]]),
                         levels = c("Control", lvl))
  # Skip if either group is too small
  n_control <- sum(dat_meta$grp == "Control")
  n_lvl <- sum(dat_meta$grp == lvl)
  if (n_control < 2 || n_lvl < 2) {
    message("  Skipping ", lvl, ": need >=2 samples in Control and group.")
    next
  }
  # Design for this pair only
  design <- model.matrix(~ 0 + grp, data = dat_meta)
  colnames(design) <- make.names(colnames(design))
  # DGEList for this pair
  y <- DGEList(counts = dat_counts)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  # Contrast: lvl - Control
  vec <- rep(0, ncol(design))
  vec[which(colnames(design) == paste0("grp", make.names(lvl)))] <- 1
  vec[which(colnames(design) == paste0("grp", make.names("Control")))] <- -1
  qlf <- glmQLFTest(fit, contrast = vec)
  tt <- topTags(qlf, n = nrow(qlf$table))$table
  out_df <- data.frame(
    Group = lvl,
    Gene = rownames(tt),
    LogFC = tt$logFC,
    LogCPM = tt$logCPM,
    PValue = tt$PValue,
    FDR = tt$FDR,
    stringsAsFactors = FALSE,
    row.names = NULL)
  out_df$score <- out_df$LogFC * (-log10(pmax(out_df$PValue, .Machine$double.xmin)))
  out_df_filtered <- out_df %>% filter(FDR < 0.05)
  out_df_filtered$Disease <- ifelse(out_df_filtered$LogFC > 0, "Up", "Down")
  out_df_filtered <- out_df_filtered %>% arrange(desc(score))
  write.csv(out_df_filtered, paste0("/COMMA/results/scrna/downstream_analysis/DGE/similarity/", lvl, "_vs_Control", ".csv"))
  dge_lists[[lvl]] <- out_df
}

# Compare Similarity of Differential Expression Lists
groups <- names(dge_lists)
cor_mat <- matrix(NA, length(groups), length(groups),
                  dimnames = list(groups, groups))
for (i in seq_along(groups)) {
  for (j in seq_along(groups)) {
    g1 <- groups[i]
    g2 <- groups[j]
    common <- intersect(dge_lists[[g1]]$Gene, dge_lists[[g2]]$Gene)
    v1 <- dge_lists[[g1]]$score[match(common, dge_lists[[g1]]$Gene)]
    v2 <- dge_lists[[g2]]$score[match(common, dge_lists[[g2]]$Gene)]
    cor_mat[i, j] <- cor(v1, v2, method = "pearson", use = "complete.obs")}}

# Plot Similarity Matrix
cor_df <- as.data.frame(as.table(cor_mat)) %>%
  setNames(c("Group1", "Group2", "Correlation")) %>%
  mutate(Group1 = factor(Group1, levels = groups),
         Group2 = factor(Group2, levels = groups))
p <- ggplot(cor_df, aes(x = Group1, y = Group2, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(    low = "#5E3C99",
    mid = "#F7F7F7",
    high = "#E66101",
                       midpoint = 0, limit = c(-1, 1),
                       name = "PearsonCorrelation") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +coord_fixed()
print(p)
ggsave("similarity_heatmap.pdf", p, width = 8, height = 7)