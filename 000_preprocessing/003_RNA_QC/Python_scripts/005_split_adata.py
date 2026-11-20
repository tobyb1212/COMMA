import pandas as pd
import re
import anndata
import scanpy as sc
import os
import random

random.seed(1234)

adata = sc.read_h5ad("/COMMA/results/scrna/count_matrix/rna_count_matrix_metadata.h5ad")

# Directory to save split .h5ad objects
output_dir = "/COMMA/results/scrna/count_matrix/study_matrices"
os.makedirs(output_dir, exist_ok=True)

# Print Metadata to .csv File
obs_csv_path = "/COMMA/results/scrna/count_matrix/study_matrices/adata_obs.csv"
adata.obs.to_csv(obs_csv_path)

# Iterate over unique studies
for study in adata.obs["Study"].dropna().unique():

    # Remove any blank spaces for filename
    study_str = str(study)
    sanitized_name = re.sub(r"\s+", "_", study_str)
    filename = f"{sanitized_name}_count_matrix.h5ad"
    filepath = os.path.join(output_dir, filename)
    
    # Subset the AnnData object by study and write to file
    subset = adata[adata.obs["Study"] == study].copy()
    subset.write(filepath)

    print(subset.shape)
    print(f"Saved {study} AnnData file to: {output_dir}")

print("Complete")
