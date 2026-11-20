library(tidyverse)

setwd("/COMMA/results/scrna/downstream_analysis/DGE")

# Files
wf <- "wilcoxon/outputs/vsmc/whole_vsmc_Cardiac_Disease_filtered.csv"
pf <- "pseudobulk/outputs/vsmc/vsmc__Cardiac_Disease_byPatient_filtered.csv"

wdf <- read.csv(wf)
pdf <- read.csv(pf)

# Merge on Group and Gene
out_df <- wdf %>%
  transmute(
    Group = group,
    Gene = names,
    Wilcoxon_Padj = pvals_adj,
    Wilcoxon_LogFC = logfoldchanges
  ) %>%
  inner_join(
    pdf %>%
      transmute(
        Group,
        Gene,
        Pseudobulk_FDR = FDR,
        Pseudobulk_LogFC = LogFC
      ),
    by = c("Group", "Gene")
  ) %>%
  filter(Wilcoxon_Padj < 0.05)

# Table 1
out_path <- "consensus/vsmc/whole_vsmc_Cardiac_Disease_filtered.csv"
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.csv(out_df, out_path, row.names = FALSE, na = "")

cat("Total consensus genes:", nrow(out_df), "\n")