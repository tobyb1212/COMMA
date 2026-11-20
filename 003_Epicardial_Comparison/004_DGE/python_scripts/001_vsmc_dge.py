import scanpy as sc 
import os 
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from collections import defaultdict
import numpy as np
import random
import pandas as pd; pd.Series.nonzero = lambda self: self.to_numpy().nonzero()
from adjustText import adjust_text

os.chdir("/COMMA/results/scrna/downstream_analysis/epicardial_coronary/DGE")
adata = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/epicardial_coronary/analysis/epicardial_integrated_final.h5ad")
adata_vsmc = adata[(adata.obs["final_vascular_cell_type"].isin(["Vascular Smooth Muscle Cell", "Atria-enriched Vascular Smooth Muscle Cell"])) & adata.obs["sample_type"].str.contains("cells", case=False, na=False)].copy()
adata = adata_vsmc.copy()

# Differential Gene Expression Analysis
sc.tl.rank_genes_groups(adata, "batch", method='wilcoxon')
sc.tl.filter_rank_genes_groups(adata, min_fold_change=1, min_in_group_fraction = 0, max_out_group_fraction = 1)
sc.pl.rank_genes_groups_dotplot(adata, key="rank_genes_groups_filtered", n_genes=10, save=f"vsmc_origin_dge.png", dendrogram = False)
# Export to .csv file
marker_genes = sc.get.rank_genes_groups_df(adata, group=None, key='rank_genes_groups_filtered')
marker_genes.to_csv(f"cell_vsmc_origin_dge.csv")

# Create Pseudobulking Dataset
def pseudobulk_with_aggregate(adata_sub, by, out_prefix, extra_meta_cols, out_dir, layer='counts'):
    """
    adata_sub: AnnData subset
    by: column in obs to aggregate by (Patient or Sample)
    out_prefix: filename prefix for saving
    extra_meta_cols: columns to extract per group for metadata
    layer: which layer to pseudobulk (None = adata.X, "counts"=raw counts)
    """
    os.makedirs(out_dir, exist_ok=True)
    agg = sc.get.aggregate(adata_sub,by=by,func="sum",layer=layer) 
    agg.X = agg.layers['sum']
    counts = agg.to_df().T  
    counts.index.name = "gene"
    counts = counts.round().astype(int)
    meta = adata_sub.obs.groupby(by).agg({c: lambda x: x.mode()[0] if len(x.mode()) > 0 else np.nan for c in extra_meta_cols})
    meta["n_cells"] = adata_sub.obs.groupby(by).size()
    counts.to_csv(f"{out_dir}/{out_prefix}_counts.csv.gz", compression="gzip")
    meta.to_csv(f"{out_dir}/{out_prefix}_metadata.csv")
    print(f"Saved pseudobulk counts → {out_dir}/{out_prefix}_counts.csv.gz")
    print(f"Saved metadata         → {out_dir}/{out_prefix}_metadata.csv")

meta_cols = ["Study", "Patient", "Sample Name", "batch"]
save_dir = "/COMMA/results/scrna/downstream_analysis/epicardial_coronary/DGE/pseudobulk"
pseudobulk_with_aggregate(adata,by="Patient",out_prefix="cell_origin",extra_meta_cols=meta_cols,out_dir=save_dir)

