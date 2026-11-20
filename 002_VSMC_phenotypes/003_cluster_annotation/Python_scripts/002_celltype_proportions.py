import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import os

os.chdir("/COMMA/results/scrna/downstream_analysis/subclustering/vsmc")
adata = sc.read_h5ad("vsmc_annotated_final.h5ad") 

# Figure 2D
ct = pd.crosstab(adata.obs['vsmc_cell_type'], adata.obs["Cardiac Disease"], normalize="index")
ct = ct[["No", "Yes"]]
colors = ["#9ecae1", "#800000"]
ax = ct.plot(kind="bar", stacked=True, figsize=(10, 5), color=colors)
ax.set_ylabel("Proportion within cluster")
ax.set_xlabel("Cluster")
ax.legend(title=f"Cardiac Disease", bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.)
plt.tight_layout()
plt.savefig("vsmc_cell_type_CD_proportions_by_cluster.pdf")
plt.close()