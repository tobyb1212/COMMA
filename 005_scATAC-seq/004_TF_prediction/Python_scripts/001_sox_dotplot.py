import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd

rna = sc.read_h5ad("//COMMA/results/scrna/downstream_analysis/subclustering/vascular/vascular_annotated_final.h5ad")

order = ["Pericyte 1", "Pericyte 2", "Pericyte 3", "Pericyte 4",
         "Atria-enriched Vascular Smooth Muscle Cell", "Vascular Smooth Muscle Cell",
         "Capillary Endothelial Cell 1", "Capillary Endothelial Cell 2", "Capillary Endothelial Cell 3",
         "Activated/Angiogenic Endothelial Cell", "Arterial Endothelial Cell 1",
         "Arterial Endothelial Cell 2", "Venous Endothelial Cell"]
        
rna.obs["vascular_cell_type"] = pd.Categorical(
    rna.obs["vascular_cell_type"],
    categories=order,
    ordered=True
)

sox_factors = [
    "SOX2",
    "SOX3",
    "SOX4",
    "SOX5",
    "SOX6",
    "SOX7",
    "SOX8",
    "SOX9",
    "SOX13",
    "SOX15",
    "SOX17",
    "SOX30",
]

# Supplementary Figure 12C 
sc.pl.dotplot(rna,var_names=sox_factors,groupby="vascular_cell_type", show=False)
plt.savefig("SOX_factors_dotplot.pdf", bbox_inches="tight")
plt.close()