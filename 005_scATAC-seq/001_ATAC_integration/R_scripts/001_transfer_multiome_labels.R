library(ArchR)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(anndata)
library(Seurat, lib.loc = "R/4.3/skylake/seurat_v5.3.0")
Sys.setenv(RETICULATE_PYTHON = "/conda/skylake/envs/r-reticulate/bin/python")
library(reticulate)
library(optparse)
library(schard)

proj <- loadArchRProject(path = '/COMMA/results/scatac/ArchR/filtered', showLogo = FALSE)
print("Loaded ArchR Project")

rna_object <- "/COMMA/results/scrna/downstream_analysis/integrated/scANVI_integrated_final_small.h5ad"
seurat_object <- schard::h5ad2seurat(rna_object)

# Filter Metadata for Only Multiome Cells
metadata <- seurat_object[[]] %>%
  filter(sequencing_method == "10X Multiome")

# Read in Barcode WhiteLists (Helps with Mapping Multiome Data)
rna_barcodes <- read_lines(gzfile("software/cellranger-arc-2.0.2/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz"))
atac_barcodes <- read_lines(gzfile("software/cellranger-arc-2.0.2/lib/python/atac/barcodes/737K-arc-v1.txt.gz"))
barcode_map <- data.frame(rna_barcode = rna_barcodes,atac_barcode = atac_barcodes,index = seq_along(rna_barcodes))

# Extract Barcodes and Sample from RNA Cell Name
rna_cells <- data.frame(cell_name = rownames(metadata)) %>%
  mutate(
    rna_barcode = str_extract(cell_name, "[^_]+$"),
    fastq_name = str_remove(cell_name, "_[^_]+$")
  ) %>%
  left_join(barcode_map, by = "rna_barcode")

# Map RNA fastq_name to ArchR linked_RNA
archr_mapping <- proj@cellColData %>%
  as.data.frame() %>%
  rownames_to_column(var = "ATAC_cell") %>%
  mutate(ATAC_sample = sub("#.*", "", ATAC_cell)) %>%
  select(ATAC_sample, linked_RNA) %>%
  group_by(linked_RNA) %>%
  summarise(atac_sample_name = first(ATAC_sample), .groups = "drop")
rna_cells <- rna_cells %>%
    left_join(archr_mapping, by = c("fastq_name" = "linked_RNA"))

# Construct New RNA Names to Match ATAC Cell Names
rna_cells <- rna_cells %>%
  filter(!is.na(atac_sample_name)) %>%
  mutate(new_cell_name = paste0(atac_sample_name, "#", rna_barcode, "-1"))
final_rna_metadata <- seurat_object[[]] %>% rownames_to_column("cell_name") %>% filter(sequencing_method == "10X Multiome")
new_cell_names <- merge(final_rna_metadata, rna_cells, by="cell_name")

# Merge ArchR Metadata with RNA Cell Name Data
archr_df <- proj@cellColData %>%
  as.data.frame() %>%
  rownames_to_column(var = "cell_name")
merged_metadata <- merge(archr_df, new_cell_names, by.x = "cell_name", by.y = "new_cell_name", all.x=T)
merged_metadata <- merged_metadata[match(rownames(proj@cellColData), merged_metadata$cell_name), ]

# Add Multiome RNA Cell Type to Metadata
proj@cellColData$multiome_rna_celltype <- NULL
proj@cellColData$multiome_rna_celltype <- merged_metadata$broad_cell_type

saveArchRProject(ArchRProj=proj,'/COMMA/results/scatac/ArchR/filtered', load = FALSE)