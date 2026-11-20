library(ArchR)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(anndata)
library(Seurat, lib.loc = "R/4.3/skylake/seurat_v5.3.0")
library(reticulate)
library(optparse)
library(schard)
library(data.table)

proj_vascular <- loadArchRProject(path = '/COMMA/results/scatac/ArchR/vascular_clustered', showLogo = FALSE)
print("Loaded ArchR Project")

library(data.table)

base_dir <- "/COMMA/data/meta_analysis/processed"

study_config <- data.table(
  Study = c("Hocker 2021", "Read 2024", "Amrute 2024", "Kanemaru 2023"),
  directory = file.path(base_dir, c("atac_hocker", "atac_hubmap", "multiome_amrute", "multiome_kanemaru")),
  technology = c("atac", "atac", "arc", "arc"))
sample_table <- merge(
  unique(data.table(Study = as.character(proj_vascular$Study),
                    Sample = as.character(proj_vascular$Sample))),
  study_config, by = "Study", all.x = TRUE)

find_fragment_file <- function(directory, sample, technology) {
  filenames <- if (technology == "arc")
    c("atac_fragments.tsv.gz", "fragments.tsv.gz")
  else
    "fragments.tsv.gz"
  candidates <- c(
    file.path(directory, sample, "outs", filenames),
    file.path(directory, sample, filenames))
  found <- candidates[file.exists(candidates)]
  if (length(found)) found[1] else NA_character_}

sample_table[, path := mapply(find_fragment_file, directory, Sample, technology)]

library(data.table)

outdir <- "/COMMA/results/scatac/multiomic_integration/scE2G/vascular/fragments"

meta <- data.table(
  Study = as.character(proj_vascular$Study),
  Sample = as.character(proj_vascular$Sample),
  cell_type = as.character(proj_vascular$vascular_cell_type),
  barcode = sub("^.*#", "", getCellNames(proj_vascular))
)

system("module load HTSlib")
system("module load SAMtools")

for (ct in unique(meta$cell_type)) {
  message("Processing: ", ct)
  tmp <- c()
  for (i in seq_len(nrow(sample_table))) {
    s <- sample_table[i]
    bc <- meta[Study == s$Study & Sample == s$Sample & cell_type == ct, barcode]
    if (!length(bc) || is.na(s$path)) next
    bcf <- tempfile(); fwrite(data.table(bc), bcf, col.names = FALSE)
    tf <- tempfile()
    prefix <- paste(gsub(" ", "_", s$Study), s$Sample, sep = "#")
    cmd <- sprintf(
      "zcat %s | awk -v p=%s 'NR==FNR{a[$1]=1;next} $4 in a{$4=p\"#\"$4;print}' OFS='\\t' %s - > %s",
      shQuote(s$path), shQuote(prefix), shQuote(bcf), shQuote(tf)
    )
    system(cmd); unlink(bcf)
    if (file.size(tf) > 0) tmp <- c(tmp, tf) else unlink(tf)
  }
  if (!length(tmp)) next
  out <- file.path(outdir, paste0(gsub("[^A-Za-z0-9_-]", "_", ct), ".fragments.tsv"))
  system(sprintf("cat %s | LC_ALL=C sort -k1,1 -k2,2n > %s",
                 paste(shQuote(tmp), collapse = " "), shQuote(out)))
  unlink(tmp)

  system2("bgzip", c("-f", out))
  system2("tabix", c("-f", "-p", "bed", paste0(out, ".gz")))
}
