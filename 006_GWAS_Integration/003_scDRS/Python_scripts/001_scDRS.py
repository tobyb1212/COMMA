import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
import matplotlib as mpl
import os

mpl.rcParams.update({
    "font.family": "Arial",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

os.chdir("/COMMA/results/scrna/downstream_analysis/subclustering/vascular/scDRS/cell_types")

order = ["Pericyte 1", "Pericyte 2", "Pericyte 3", "Pericyte 4",
         "Atria-enriched Vascular Smooth Muscle Cell", "Vascular Smooth Muscle Cell",
         "Capillary Endothelial Cell 1", "Capillary Endothelial Cell 2",
         "Capillary Endothelial Cell 3", "Activated/Angiogenic Endothelial Cell",
         "Arterial Endothelial Cell 1", "Arterial Endothelial Cell 2",
         "Venous Endothelial Cell"]

df = pd.read_csv("CAD.scdrs_group.vascular_cell_type", sep="\t")

df = df.set_index("group").loc[order].reset_index()
df["logP"] = -np.log10(df["assoc_mcp"])
df["sig"] = df["assoc_mcp"] < 0.05 / len(df)
df["y"] = np.arange(len(df))

lim = abs(df["assoc_mcz"]).max()
norm = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
cmap = plt.get_cmap("RdBu_r")
size = 50 + 100 * df["logP"]

fig, ax = plt.subplots(figsize=(8, 5.5))

for sig, d in df.groupby("sig"):
    ax.scatter(d["assoc_mcz"], d["y"], s=size[d.index],
               c=d["assoc_mcz"] if sig else "none",
               edgecolors=cmap(norm(d["assoc_mcz"])),
               cmap=cmap if sig else None, norm=norm if sig else None,
               linewidth=1.2)

ax.set_yticks(df["y"], df["group"])
ax.invert_yaxis()  # Pericyte 1 at top

ax.axvline(0, ls="--", lw=1, alpha=.6)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
fig.colorbar(sm, ax=ax, pad=.02, label="scDRS association Z-score")

ax.legend(handles=[
    Line2D([0],[0], marker="o", ls="", mfc="black", mec="black",
           label="Bonferroni significant"),
    Line2D([0],[0], marker="o", ls="", mfc="none", mec="black",
           label="Not significant")
], frameon=False)

ax.set(xlabel="CAD scDRS association Z-score", ylabel="",
       title="CAD genetic association across vascular cell types")
ax.spines[["top","right"]].set_visible(False)

# Figure 6C
plt.tight_layout()
plt.savefig("CAD_scDRS_celltype_dotplot.png", dpi=300, bbox_inches="tight")
plt.close()