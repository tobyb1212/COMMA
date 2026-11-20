library(data.table)
library(edgeR)
library(tidyverse)
library(ggrepel)

# Define Directory & Conditions to Compare
base_dir <- "/COMMA/results/scrna/downstream_analysis/epicardial_coronary/DGE/pseudobulk"
counts_fp <- file.path(base_dir, "cell_origin_counts.csv.gz")
metadata_fp <- file.path(base_dir, "cell_origin_metadata.csv")
group_col <- "batch"

# Read In Counts
counts_dt <- fread(counts_fp, data.table = FALSE)
rownames(counts_dt) <- counts_dt[[1]]
counts_dt[[1]] <- NULL
counts_mat <- as.matrix(counts_dt)

# Align Metadata to Counts
meta <- read.csv(metadata_fp, row.names = 1, check.names = FALSE)
meta <- meta[colnames(counts_mat), , drop = FALSE]
meta$grp <- factor(meta[[group_col]])

# Prepare edgeR Design Matrix
design <- model.matrix(~ 0 + grp, data = meta)
colnames(design) <- make.names(colnames(design))

# Perform edgeR Preparation
y <- DGEList(counts = counts_mat)
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# Perform edgeR DGE Over Each Group (one-vs-rest)
all_contrasts <- list()

for (j in seq_len(ncol(design))) {

  level_name <- levels(meta$grp)[j]
  n_in <- sum(meta$grp == level_name)
  n_others <- nrow(meta) - n_in
  if (n_in < 2 || n_others < 2) next

  contrast <- rep(-1 / (ncol(design) - 1), ncol(design))
  contrast[j] <- 1

  qlf <- glmQLFTest(fit, contrast = contrast)
  tt <- topTags(qlf, n = Inf)$table

  out_df <- data.frame(
    Group = level_name,
    Gene = rownames(tt),
    LogFC = tt$logFC,
    LogCPM = tt$logCPM,
    FDR = tt$FDR
  )

  # Apply filtering: magnitude of fold change > 0.5 and FDR < 0.05
  filtered <- out_df %>% filter(LogFC > 0.5, FDR < 0.05)
  all_contrasts[[level_name]] <- filtered
}

# Combine Output Files and Write
combined <- bind_rows(all_contrasts)
out_path <- file.path(base_dir, "outputs", "cell_origin_filtered.csv")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.csv(combined, out_path, row.names = FALSE, quote = FALSE)

cat("Total filtered genes:", nrow(combined), "\n")

# Volcano Plot
all_genes <- topTags(qlf, n = Inf)$table %>%
  rownames_to_column("gene")

volc <- all_genes %>%
  transmute(
    EFFECTSIZE = logFC,
    p.adj = -log10(pmax(FDR, .Machine$double.xmin)),
    gene,
    Label = gene %in% c("MGP", "LUM", "VIM", "ADAMTS9", "PDE1C")
  )

# Figure 3D
ggplot(volc, aes(EFFECTSIZE, p.adj)) +
  geom_point(color = "#878787", alpha = 0.2) +
  geom_point(data = filter(volc, p.adj > -log10(0.05), EFFECTSIZE > 1),
             color = "#E9A23B") +
  geom_point(data = filter(volc, p.adj > -log10(0.05), EFFECTSIZE < -1),
             color = "#2A9D8F") +
  geom_vline(xintercept = 0, colour = "black", linetype = "longdash") +
  scale_x_continuous(limits = c(-12.5, 12.5), breaks = seq(-10, 10, 5)) +
  scale_y_continuous(limits = c(0, 15)) +
  labs(
    x = expression(log[2] * " Fold Change"),
    y = expression(-log[10] * " adjusted p-value")
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_text_repel(
    data = filter(volc, Label),
    aes(label = gene),
    box.padding = 0.5,
    point.padding = 0.3,
    force = 2,
    max.overlaps = Inf,
    segment.color = "black",
    segment.size = 0.4,
    segment.alpha = 1,
    min.segment.length = 0
  )

ggsave(
  filename = file.path(base_dir, "outputs", "volcano.pdf"),
  plot = last_plot(),
  height = 80, width = 100, units = "mm"
)