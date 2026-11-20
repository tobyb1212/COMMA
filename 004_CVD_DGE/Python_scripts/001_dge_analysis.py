import scanpy as sc
import pandas as pd
import os
import math
from collections import defaultdict
import numpy as np
import random

random.seed(1234)

vsmc_adata = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/subclustering/vsmc/vsmc_annotated_final.h5ad")

os.chdir("/COMMA/results/scrna/downstream_analysis/DGE/wilcoxon/outputs/vsmc")

for condition in ['Cardiac Disease']:
    condition_safe = condition.replace(" ", "_")
    vsmc_adata_condition = vsmc_adata[~vsmc_adata.obs[condition].isna(), :].copy()
    # Differential Gene Expression for Each Condition
    sc.tl.rank_genes_groups(vsmc_adata_condition, condition, method='wilcoxon')
    sc.tl.filter_rank_genes_groups(vsmc_adata_condition, min_fold_change=0.5)
    marker_df = sc.get.rank_genes_groups_df(vsmc_adata_condition, group=None, key='rank_genes_groups_filtered')
    marker_df.to_csv(f"whole_vsmc_{condition_safe}_unfiltered_genes.csv")
    # Keep only valid studies for per-study analysis
    valid_studies = []
    for Study in vsmc_adata_condition.obs["Study"].unique():
        group_counts = vsmc_adata_condition.obs.loc[
            vsmc_adata_condition.obs["Study"] == Study, condition
        ].value_counts()
        if len(group_counts) >= 2 and (group_counts >= 3).all():
            valid_studies.append(Study)
    vsmc_adata_condition = vsmc_adata_condition[ vsmc_adata_condition.obs["Study"].isin(valid_studies)].copy()
    # Filter Genes that aren't DEGs in half of the studies when compared individually
    counts = defaultdict(lambda: defaultdict(int))
    logfcs = defaultdict(lambda: defaultdict(list))
    processed_studies = set()
    for Study in vsmc_adata_condition.obs["Study"].unique().tolist():
        try:
            adata_subset = vsmc_adata_condition[vsmc_adata_condition.obs["Study"] == Study].copy()
            adata_subset.obs[condition] = adata_subset.obs[condition].astype("category").cat.remove_unused_categories()
            sc.tl.rank_genes_groups(adata_subset, groupby=condition, method="wilcoxon")
            sc.tl.filter_rank_genes_groups(adata_subset, min_fold_change=0.5, min_in_group_fraction=0, max_out_group_fraction = 1)
            df = sc.get.rank_genes_groups_df(adata_subset, group=None, key="rank_genes_groups_filtered")
            df = df[df["pvals_adj"] < 0.05].copy()
        except Exception as e:
            print(f"Skipping Study={Study} due to error: {e}")
            continue
        if df is None or df.empty:
            continue
        processed_studies.add(Study)
        for grp, gdf in df.groupby("group"):
            overall = marker_df[marker_df["group"] == grp].set_index("names")["logfoldchanges"]
            for _, r in gdf.iterrows():
                g = r["names"]
                # Only count genes with the same fold-change direction as overall
                if pd.isna(g) or g not in overall.index or pd.isna(r["logfoldchanges"]):
                    continue
                if r["logfoldchanges"] * overall.loc[g] <= 0:
                    continue
                counts[grp][g] += 1
                logfcs[grp][g].append(r["logfoldchanges"])
    rows = []
    for cl in counts:
        for g, n in counts[cl].items():
            avg_lfc = sum(logfcs[cl][g]) / len(logfcs[cl][g])
            rows.append([cl, g, n, avg_lfc])
    # Filter Final Markers by Criteria
    out = pd.DataFrame(rows, columns=["cluster", "gene", "n_studies", "avg_logFC"])
    n_processed = len(processed_studies)
    threshold = math.ceil(n_processed / 2.0)
    out_filtered = out[out["n_studies"] >= threshold].copy()
    out_filtered = out_filtered[out_filtered["gene"].notna() & (out_filtered["gene"].astype(str).str.strip() != "")].copy()
    keep_pairs = set(zip(out_filtered["cluster"].astype(str), out_filtered["gene"].astype(str)))
    mask = marker_df.apply(lambda r: (str(r['group']), str(r['names'])) in keep_pairs, axis=1)
    filtered_df = marker_df[mask].copy()
    filtered_df.to_csv(f"whole_vsmc_{condition_safe}_filtered.csv")
    # Plot DotPlot for Top 5 Genes in Each Condition
    top_genes = (filtered_df.sort_values(["group", "scores"], ascending=[True, False]).groupby("group").head(5))
    genes_dict = {grp: df["names"].tolist() for grp, df in top_genes.groupby("group")}
    all_genes = list(dict.fromkeys(top_genes["names"].tolist()))
    if genes_dict:
        sc.pl.dotplot(vsmc_adata_condition, var_names=genes_dict, groupby=condition,
            figsize=(8, 6), save=f"_{condition_safe}_filtered_dotplot.png")

# Supplemental Figure 8A
vsmc_adata.obs["Cardiac Disease x sequencing_method"] = (vsmc_adata.obs["Cardiac Disease"].astype(str) + " x " + vsmc_adata.obs["sequencing_method"].astype(str))
sc.pl.dotplot(vsmc_adata,var_names=["ADAMTS9", "ADAMTS9-AS2"],groupby="Cardiac Disease x sequencing_method",swap_axes=True, save=f"ADAMTS9_locus_plot.pdf")