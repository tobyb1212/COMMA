library(data.table)
library(edgeR)
library(tidyverse)
library(ggrepel)

# Define Directory & Conditions to Compare
base_dir <- "/COMMA/results/scrna/downstream_analysis/DGE/pseudobulk"
fname <- "vsmc__Cardiac_Disease_byPatient"

counts_fp <- file.path(base_dir, "vsmc", paste0(fname, "_counts.csv.gz"))
metadata_fp <- file.path(base_dir, "vsmc", paste0(fname, "_metadata.csv"))

# Read In Counts
counts_dt <- fread(counts_fp, data.table = FALSE)
rownames(counts_dt) <- counts_dt[[1]]
counts_dt[[1]] <- NULL
counts_mat <- as.matrix(counts_dt)
mode(counts_mat) <- "numeric"

# Align Metadata to Counts
meta <- read.csv(metadata_fp, row.names = 1, check.names = FALSE)
stopifnot(all(colnames(counts_mat) %in% rownames(meta)))
meta <- meta[colnames(counts_mat), , drop = FALSE]
meta$grp <- factor(meta[["Cardiac Disease"]])
meta$Study <- factor(meta$Study)

# Prepare edgeR Design Matrix
design <- model.matrix(~ 0 + grp + Study, data = meta)
colnames(design) <- make.names(colnames(design))

# Perform edgeR Preparation
y <- DGEList(counts = counts_mat)
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# Perform edgeR DGE Iteratively Over Each Group
all_contrasts <- list()
group_cols <- grep("^grp", colnames(design))

for (j in seq_along(group_cols)) {

  level_name <- levels(meta$grp)[j]
  n_in <- sum(meta$grp == level_name)
  n_others <- nrow(meta) - n_in
  if (n_in < 2 || n_others < 2) next

  contrast <- rep(0, ncol(design))
  contrast[group_cols] <- -1 / (length(group_cols) - 1)
  contrast[group_cols[j]] <- 1

  qlf <- glmQLFTest(fit, contrast = contrast)
  tt <- topTags(qlf, n = Inf)$table

  out_df <- data.frame(
    Group = level_name,
    Gene = rownames(tt),
    LogFC = tt$logFC,
    LogCPM = tt$logCPM,
    FDR = tt$FDR,
    group = level_name
  )

  # Apply filtering: magnitude of fold change > 0.5 and FDR < 0.05
  all_contrasts[[level_name]] <- out_df %>%
    filter(LogFC > 0.5, FDR < 0.05)
}

# Combine Output Files and Write
combined <- bind_rows(all_contrasts)
out_path <- file.path(base_dir, "outputs/vsmc", paste0(fname, "_filtered.csv"))
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.csv(combined, out_path, row.names = FALSE, quote = FALSE)

cat("Total filtered genes:", nrow(combined), "\n")

# Figure 4A
all_genes <- topTags(qlf, n = Inf)$table %>%
  rownames_to_column("Gene")
interested_genes <- c(
  "ADAMTS9", "PTP4A3", "RGL3", "COL4A1", "BTBD9",
  "PDE10A", "ADAMTS9-AS2", "ITIH5", "TMEM51", "ZNF704"
)
volc <- all_genes %>%
  transmute(
    EFFECTSIZE = logFC,
    p.adj = -log10(pmax(FDR, .Machine$double.xmin)),
    gene = Gene,
    Label = Gene %in% interested_genes
  )
ggplot(volc, aes(EFFECTSIZE, p.adj)) +
  geom_point(color = "#878787", alpha = 0.2) +
  geom_point(
    data = filter(volc, p.adj > -log10(0.05), EFFECTSIZE > 0.5),
    color = "red"
  ) +
  geom_point(
    data = filter(volc, p.adj > -log10(0.05), EFFECTSIZE < -0.5),
    color = "blue"
  ) +
  geom_vline(xintercept = 0, colour = "black", linetype = "longdash") +
  scale_x_continuous(limits = c(-5, 5), breaks = seq(-5, 5, 1)) +
  scale_y_continuous(limits = c(0, 10)) +
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
  filename = "/COMMA/results/scrna/downstream_analysis/DGE/consensus/Vascular_Smooth_Muscle_Cell/volcano.pdf",
  plot = last_plot(),
  height = 80, width = 100, units = "mm"
)