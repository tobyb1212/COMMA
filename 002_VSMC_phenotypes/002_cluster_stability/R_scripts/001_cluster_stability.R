library(scclusteval)
library(tidyverse)
library(patchwork)
library(SummarizedExperiment)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SingleCellExperiment)
library(anndata)
library(Seurat, lib.loc = "R/4.3/skylake/seurat_v5.3.0")
Sys.setenv(RETICULATE_PYTHON = "/conda/skylake/envs/r-reticulate/bin/python")
library(harmony)
library(reticulate)
library(optparse)
library(schard)
library(clue)  
library(purrr)

setwd("/COMMA/results/scrna/downstream_analysis/subclustering/cluster_stability")

# Convert .h5ad Object to Anndata
seurat_object <- schard::h5ad2seurat("/COMMA/results/scrna/downstream_analysis/subclustering/vsmc/scVI_vsmc_integrated.h5ad")
seurat_counts <- schard::h5ad2seurat("/COMMA/results/scrna/downstream_analysis/subclustering/vsmc/scVI_vsmc_integrated_counts.h5ad")
seurat_object@assays$RNA@counts <- seurat_counts@assays$RNA@counts
saveRDS(seurat_object, "/COMMA/results/scrna/downstream_analysis/subclustering/cluster_stability.rds")

# Subsample Data & Iterations of Subsampling Done with PyflowSeurat
system("sbatch software/pyflow_seurat_parameter_adapted/pyflow-scBoot.sh")

# Read in Gathered Full Sample & Gathered Sub Sample
subsample_idents<- readRDS("software/pyflow_seurat_parameter_adapted/gather_subsample.rds")
fullsample_idents<- readRDS("software/pyflow_seurat_parameter_adapted/gather_full_sample.rds")

subsample_idents <- subsample_idents %>%
  mutate(
    original_ident = map(original_ident, ~ {
      names(.x) <- paste0("cell_", seq_along(.x))
      .x
    }),
    recluster_ident = map2(recluster_ident, original_ident, ~ {
      names(.x) <- names(.y)
      .x
    })
  )

subsample_idents_list <- subsample_idents %>% 
  group_by(pc, resolution, k.param) %>% 
  nest()

# Assign Stable Clusters
stable_clusters<- subsample_idents_list %>%
  mutate(stable_cluster = map(data, ~ AssignStableCluster(.x$original_ident,
                                                          .x$recluster_ident,
                                                          jaccard_cutoff = 0.6,
                                                          method = "jaccard_percent", 
                                                          percent_cutoff = 0.6)))

# Parameter Scatter Plot
# Supplementary Figure 6A
ParameterSetScatterPlot(stable_clusters = stable_clusters,
                        fullsample_idents = fullsample_idents,
                        x_var = "k.param",
                        y_var = "number",
                        facet_rows = "resolution",
                        facet_cols = "pc") + theme_bw()
ggsave(plot = last_plot(), filename = "parameter_scatter.pdf", width = 10, height = 10)

# Plot Jacard RainCloud Plot
# Supplementary FIgure 6B
JaccardRainCloudPlot(subsample_idents_list$data[[54]]$original_ident,
                          subsample_idents_list$data[[54]]$recluster_ident) + 
        geom_hline(yintercept = c(0.6, 0.75), linetype = 2) +
        xlab("Cluster ID (Neighbours=45 Resolution=0.3 Components=30)") 
ggsave(plot = last_plot(), filename = "parameter_raincloud_k45_res0.3_comp30.pdf", height = 6, width = 10)
