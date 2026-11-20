# Rscript
library(tidyverse)

base <- "/COMMA/results/scrna/downstream_analysis/epicardial_coronary/DGE"

# Read files
wdf <- read.csv(
  file.path(base, "cell_vsmc_origin_dge.csv"),
  check.names = FALSE
)

wdf <- wdf[, names(wdf) != "", drop = FALSE]

pdf <- read.csv(
  file.path(base, "pseudobulk/outputs/cell_origin_filtered.csv"),
  check.names = FALSE
)

# Prepare Wilcoxon
wdf2 <- wdf %>%
  transmute(
    Group = group,
    Gene = names,
    Wilcoxon_Padj = as.numeric(pvals_adj),
    Wilcoxon_LogFC = as.numeric(logfoldchanges),
    Group_u = toupper(str_squish(group)),
    Gene_u = toupper(str_trim(names))
  )

# Prepare pseudobulk
pdf2 <- pdf %>%
  transmute(
    Pseudobulk_LogFDR = as.numeric(FDR),
    Pseudobulk_LogFC = as.numeric(LogFC),
    Group_u = toupper(str_squish(Group)),
    Gene_u = toupper(str_trim(Gene))
  )

# Check group labels
cat("Wilcoxon groups:\n")
print(unique(wdf2$Group_u))

cat("\nPseudobulk groups:\n")
print(unique(pdf2$Group_u))

# Merge on group and gene
consensus <- inner_join(
  wdf2,
  pdf2,
  by = c("Group_u", "Gene_u")
) %>%
  select(
    Group,
    Gene,
    Wilcoxon_Padj,
    Wilcoxon_LogFC,
    Pseudobulk_LogFDR,
    Pseudobulk_LogFC
  )

consensus <- consensus %>% filter(!is.na(Wilcoxon_Padj) & Wilcoxon_Padj < 0.05)

# Supplemental Table 8
out_path <- file.path(base, "cell_consensus_genes.csv")

write.csv(consensus, out_path, row.names = FALSE, na = "")

cat("\nConsensus genes:", nrow(consensus), "\n")
cat("Saved to:", out_path, "\n")