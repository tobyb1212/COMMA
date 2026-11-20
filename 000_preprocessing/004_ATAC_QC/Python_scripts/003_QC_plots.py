import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np

atac = sc.read_h5ad("/COMMA/results/scatac/atac_qc.h5ad")

# Supplementary Figure 10
plt.figure(figsize=(6, 3.2))
ax = sc.pl.violin(atac,["TSSEnrichment"], groupby='Study', stripplot=False, show=False, rotation=90)
ax.axhline(y=3, color='red', linestyle='dotted', linewidth=1)
plt.ylim(0, 20) 
plt.savefig("cell_metrics/TSS_enrichment.png", dpi=300, bbox_inches="tight")
plt.close()
plt.figure(figsize=(6, 3.2))
ax = sc.pl.violin(atac,["DoubletScore"], groupby='Study', stripplot=False, show=False, rotation=90)
plt.ylim(0, 20) 
plt.savefig("cell_metrics/doublet_violin.png", dpi=300, bbox_inches="tight")
plt.close()
plt.figure(figsize=(6, 3.2))
atac.obs["logFrags"] = np.log10(atac.obs["nFrags"])
ax = sc.pl.violin(atac,["logFrags"], groupby='Study', stripplot=False, show=False, rotation=90)
plt.ylim(0, 7)
ax.axhline(y=3, color='red', linestyle='dotted', linewidth=1)
plt.savefig("cell_metrics/logFrags.png", dpi=300, bbox_inches="tight")
plt.close()