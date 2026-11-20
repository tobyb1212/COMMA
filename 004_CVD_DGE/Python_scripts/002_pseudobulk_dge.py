import scanpy as sc
import pandas as pd
import os
import math
from collections import defaultdict
import numpy as np
import random

random.seed(1234)

# Function to Pseudobulk & Save Aggregated Counts with Metadata
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

# Perform VSMC Pseudobulking Analysis
conditions = ['Cardiac Disease', 'Primary Cardiac Disease']
vsmc_adata = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/subclustering/vsmc/vsmc_annotated_final.h5ad")

# Perform Pseudobulking for Whole VSMCs
for cond in conditions:
    save_dir = "/COMMA/results/scrna/downstream_analysis/DGE/pseudobulk/vsmc"
    vsmc_adata_condition = vsmc_adata[~vsmc_adata.obs[cond].isna(), :].copy()
    # Rule: Region → by Sample, else by Patient
    if cond in ["Region", "Atria_Ventricle"]:
            group_by = "Sample Name"
            group_label = "bySample"
    else:
            group_by = "Patient"
            group_label = "byPatient"
    out_prefix = f"vsmc__{cond.replace(' ', '_')}_{group_label}"
    vsmc_adata_condition = vsmc_adata_condition[vsmc_adata_condition.obs[group_by].map(vsmc_adata_condition.obs[group_by].value_counts()) >= 30].copy()
    meta_cols = ["Study", "Patient", "Sample Name", cond]
    pseudobulk_with_aggregate(vsmc_adata_condition,by=group_by,out_prefix=out_prefix,extra_meta_cols=meta_cols,out_dir=save_dir)

