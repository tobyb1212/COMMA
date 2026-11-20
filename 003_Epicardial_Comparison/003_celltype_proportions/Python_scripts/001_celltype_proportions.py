import sccoda
import scanpy as sc 
import os 
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sccoda.util import cell_composition_data
from sccoda.util import comp_ana as mod
from sccoda.util import cell_composition_data as dat
from sccoda.util import data_visualization as viz

os.chdir("/COMMA/results/scrna/downstream_analysis/epicardial_coronary/analysis")
adata = sc.read_h5ad("no_uns_epicardial_integrated_final_small.h5ad")

# Define Variables
sample_col = "Sample Name"  
cell_type_col = "final_vascular_cell_type"
condition_col = "batch"  
adata = adata[adata.obs["final_vascular_cell_type"] != "Unknown"].copy()
adata.obs["final_vascular_cell_type"] = adata.obs["final_vascular_cell_type"].cat.remove_unused_categories()
adata.obs["Sample Name"] = adata.obs["Sample Name"].astype(str)
cov_df = (adata.obs[[sample_col, "batch"]].drop_duplicates()
    .set_index(sample_col).sort_index())
data_scanpy = dat.from_scanpy(adata, cell_type_identifier=cell_type_col,
    sample_identifier=sample_col, covariate_df=cov_df)

# Data from 
viz.boxplots(data_scanpy,feature_name="batch", plot_facets=False, y_scale="relative", add_dots=False)
plt.savefig("test_proportions.png")

# Construct Model
model_all = mod.CompositionalAnalysis(data_scanpy, formula=condition_col, reference_cell_type="automatic")
all_results = model_all.sample_hmc()
all_results.summary()

# Plot Only Significant Findings
cred = all_results.credible_effects(est_fdr=0.01)
cred_df = cred.reset_index()
credible_celltypes = (cred_df.loc[cred_df["Final Parameter"], "Cell Type"].unique().tolist())

def add_credible_brackets(ax, df, x_col, y_col, credible_labels, star="*", fontsize=14, bar_width=0.35, tick_height_linear=0.03, tick_height_log=0.08, offset_linear=0.07, offset_log=0.20):
    xticks = ax.get_xticks()
    xticklabels = [t.get_text() for t in ax.get_xticklabels()]
    is_log = ax.get_yscale() == "log"
    for x, label in zip(xticks, xticklabels):
        if label in credible_labels:
            ymax = df.loc[df[x_col] == label, y_col].max()
            if is_log:
                y = ymax * (1 + offset_log)
                tick = y * tick_height_log
            else:
                y = ymax * (1 + offset_linear)
                tick = tick_height_linear * ymax
            ax.plot([x - bar_width / 2, x + bar_width / 2], [y, y], lw=1.8, c="black", clip_on=False)
            ax.plot([x - bar_width / 2, x - bar_width / 2], [y, y - tick], lw=1.8, c="black", clip_on=False)
            ax.plot([x + bar_width / 2, x + bar_width / 2], [y, y - tick], lw=1.8, c="black", clip_on=False)
            ax.text(x, y + tick * 0.6, star, ha="center", va="bottom", fontsize=fontsize, fontweight="bold", color="black", clip_on=False)

def compute_proportions(adata):
    obs = adata.obs.copy()
    counts = obs.groupby(["Sample Name", "final_vascular_cell_type"]).size().reset_index(name="n_cells")
    counts = counts.merge(obs[["Sample Name", "batch"]].drop_duplicates(), on="Sample Name", how="left")
    totals = obs.groupby("Sample Name").size().reset_index(name="total_cells")
    df = counts.merge(totals, on="Sample Name")
    df["proportion"] = df["n_cells"] / df["total_cells"]
    return df

df_all = compute_proportions(adata)
df_all = df_all.copy()
df_all["proportion_log"] = df_all["proportion"] + 1e-6
order = ["Pericyte 1", "Pericyte 2", "Pericyte 3", "Pericyte 4", "Atria-enriched Vascular Smooth Muscle Cell", "Vascular Smooth Muscle Cell", "Capillary Endothelial Cell 1", "Capillary Endothelial Cell 2", "Capillary Endothelial Cell 3", "Activated/Angiogenic Endothelial Cell", "Arterial Endothelial Cell 1", "Arterial Endothelial Cell 2", "Venous Endothelial Cell"]
df_all["final_vascular_cell_type"] = pd.Categorical(df_all["final_vascular_cell_type"], categories=order, ordered=True)
df_all = df_all.sort_values("final_vascular_cell_type")

# Figure 3C
fig, ax = plt.subplots(figsize=(14, 6))
sns.boxplot(data=df_all, x="final_vascular_cell_type", y="proportion_log", hue="batch", showfliers=False, ax=ax)
sns.stripplot(data=df_all, x="final_vascular_cell_type", y="proportion_log", hue="batch", dodge=True, color="black", alpha=0.7, size=1, linewidth=0.5, edgecolor="black", ax=ax)
ax.set_yscale("log")
ax.set_ylabel("Proportion (log scale)")
ax.set_xlabel("Epicardial-Coronary Cell Type")
ax.set_title("Cell Type Proportions (Log scale)")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:2], labels[:2], title="Origin")
add_credible_brackets(ax=ax, df=df_all, x_col="final_vascular_cell_type", y_col="proportion_log", credible_labels=credible_celltypes)
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("significant_celltype_proportions_log.pdf", dpi=300, bbox_inches="tight")
plt.show()

