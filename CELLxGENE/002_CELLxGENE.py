import pandas as pd
import numpy as np
import scanpy as sc
import os
import argparse
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
import re
import gzip
from pathlib import Path
import shutil

vascular_adata = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/subclustering/vascular/vascular_annotated_final.h5ad")
vsmc_adata = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/subclustering/vsmc/vsmc_annotated_final.h5ad")
atac_adata = sc.read_h5ad("/COMMA/results/data_submission/raw/genescore_matrix.h5ad")

# Filter out dbGaP Samples
dbgap_studies = ["Tucker 2020", "Hocker 2021", "Chaffin 2022", "Simonson 2023", "Read 2024"]
vascular_adata = vascular_adata[~vascular_adata.obs["Study"].isin(dbgap_studies)].copy()
vsmc_adata = vsmc_adata[~vsmc_adata.obs["Study"].isin(dbgap_studies)].copy()
atac_adata = atac_adata[~atac_adata.obs["Study"].isin(dbgap_studies)].copy()

# Convert ENSEMBL IDs to Match v114
gtf46 = "/COMMA/reference_genome/GENCODE/GRCh38.p14.GENCODE_release46/gencode.v46.basic.annotation_filtered.gtf"
gtf48 = "/COMMA/reference_genome/GENCODE/CZI_CELLxGENE/gencode.v48.primary_assembly.annotation.gtf.gz"
genes = {}
with open(gtf46) as f:
    for line in f:
        if line.startswith("#") or line.split("\t")[2] != "gene": continue
        attrs = dict(re.findall(r'(\w+) "([^"]+)"', line.split("\t")[8]))
        genes.setdefault(attrs["gene_name"], set()).add(attrs["gene_id"].split(".")[0])

gene_map = {k: next(iter(v)) for k, v in genes.items() if len(v) == 1}
with gzip.open(gtf48, "rt") as f:
    valid_ids = {re.search(r'gene_id "([^"]+)"', line).group(1).split(".")[0] for line in f if not line.startswith("#") and "\tgene\t" in line}

# Convert to GENCODE v46 ENSEMBL IDs and confirm v48 ENSEMBL ID Match. Drop IDs not in ENSEMBL ID Match 
for name, adata in [("Vascular", vascular_adata), ("VSMC", vsmc_adata)]:
    ids = pd.Index([x.split(".")[0] if x.startswith("ENSG") else gene_map.get(x, x) for x in adata.var_names])
    keep = ids.isin(valid_ids) & ~ids.duplicated(keep=False)
    print(f"{name}: keeping {keep.sum()}/{len(keep)} genes, dropping {(~keep).sum()}")
    adata._inplace_subset_var(keep)
    adata.var_names = ids[keep]

# Convery ATAC Genes to Ensembl IDs
genes48 = {}
with gzip.open(gtf48, "rt") as f:
    for line in f:
        if line.startswith("#") or line.split("\t")[2] != "gene": continue
        attrs = dict(re.findall(r'(\w+) "([^"]+)"', line.split("\t")[8]))
        genes48.setdefault(attrs["gene_name"], set()).add(attrs["gene_id"].split(".")[0])

gene_map48 = {k: next(iter(v)) for k, v in genes48.items() if len(v) == 1}
ids = pd.Index([x.split(".")[0] if x.startswith("ENSG") else gene_map48.get(x, x) for x in atac_adata.var_names])
keep = ids.isin(valid_ids) & ~ids.duplicated(keep=False)
print(f"ATAC: keeping {keep.sum()}/{len(keep)} genes, dropping {(~keep).sum()}")
print("Unmapped:", atac_adata.var_names[~keep][:20].tolist())
atac_adata._inplace_subset_var(keep)
atac_adata.var_names = ids[keep]

SEX_MAP = {"Male": "PATO:0000384", "Female": "PATO:0000383"}
SUSPENSION_MAP = {"Nuclei": "nucleus", "Cells": "cell", "CD45+ cells": "cell", "CD45- cells": "cell"}
REGION_MAP = {"Apex": "UBERON:0002084", "LV": "UBERON:0002084", "LA": "UBERON:0002079", "RA": "UBERON:0002078", "Septum": "UBERON:0002094", "RV": "UBERON:0002080", "SAN": "UBERON:0002351", "AVN": "UBERON:0002352"}
DISEASE_MAP = {"Control": "PATO:0000461", "DCM": "MONDO:0005021", "IHD": "MONDO:0024644", "NICM": "MONDO:0004994", "AF": "MONDO:0004981", "AMI": "MONDO:0005068", "CAD": "MONDO:0005010", "Cardiac Arrest": "MONDO:0000745"}
RNA_ASSAY_MAP = {"10X V3": "EFO:0009922", "10X Chromium 5ʹ V1.1": "EFO:0011025", "10X V2": "EFO:0009899", "10X Multiome": "EFO:0030059"}
ATAC_ASSAY_MAP = {"10X Multiome": "EFO:0030059", "Combinatorial Indexing": "EFO:0022045"}
CELL_TYPE_MAP = {"Pericyte": "CL:0000669", "Pericyte 1": "CL:0000669", "Pericyte 2": "CL:0000669", "Pericyte 3": "CL:0000669", "Pericyte 4": "CL:0000669", "Capillary Endothelial Cell": "CL:0002144", "Capillary Endothelial Cell 1": "CL:0002144", "Capillary Endothelial Cell 2": "CL:0002144", "Capillary Endothelial Cell 3": "CL:0002144", "Arterial Endothelial Cell": "CL:1000413", "Arterial Endothelial Cell 1": "CL:1000413", "Arterial Endothelial Cell 2": "CL:1000413", "Venous Endothelial Cell": "CL:0002543", "Activated/Angiogenic Endothelial Cell": "CL:0000115", "Vascular Smooth Muscle Cell": "CL:0000359", "Atria-enriched Vascular Smooth Muscle Cell": "CL:0000359"}
ETHNICITY_MAP = {"Caucasian": "HANCESTRO:0005", "African American": "HANCESTRO:0568", "Asian": "HANCESTRO:0008", "South Asian": "HANCESTRO:0006", "Hispanic": "HANCESTRO:0014"}
adatas = [(vascular_adata, "COMMA sc/snRNA-seq", "rna"), (vsmc_adata, "sc/snRNA-seq Microvascular SMC Atlas", "rna"), (atac_adata, "COMMA scATAC-seq", "atac")]
DOI_MAP = {"Hocker 2021": "10.1126/sciadv.abf1444", "Litvinukova 2020": "10.1038/s41586-020-2797-4", "Tucker 2020": "10.1161/CIRCULATIONAHA.119.045401", "Rao 2021": "10.1007/s00395-021-00897-1", "Chaffin 2022": "10.1038/s41586-022-04817-8","Koenig 2022": "10.1038/s44161-022-00028-6", "Amrute 2023": "10.1038/s44161-023-00260-8", "Kanemaru 2023": "10.1038/s41586-023-06311-1","Simonson 2023": "10.1016/j.celrep.2023.112086", "Amrute 2024": "10.1038/s41586-024-08008-5","Hill 2024": "10.1038/s41467-024-54296-w","Linna-Kuosmanen 2024": "10.1016/j.xcrm.2024.101556", "Read 2024": "10.1038/s42003-024-06582-y",}

for adata, title, modality in adatas:
    adata.uns['title'] = title
    adata.uns['organism_ontology_term_id'] = "NCBITaxon:9606"
    if modality == 'rna':
        adata.raw = sc.AnnData(X=adata.layers['counts'].astype(np.float32).copy(),var=adata.var.copy())
    adata.layers.clear()
    adata.obs['donor_id'] = adata.obs['Patient']
    adata.obs['suspension_type'] = adata.obs['sample_type'].map(SUSPENSION_MAP)
    adata.obs['tissue_type'] = 'tissue'
    adata.obs['sex_ontology_term_id'] = adata.obs['Sex'].map(SEX_MAP)
    adata.obs['development_stage_ontology_term_id'] = 'HsapDv:0000258'
    adata.obs['tissue_ontology_term_id'] = adata.obs['Region'].map(REGION_MAP)
    adata.obs['cell_type_ontology_term_id'] = adata.obs['vascular_cell_type'].map(CELL_TYPE_MAP) if adata is not vsmc_adata else 'CL:0000359'
    adata.obs['disease_ontology_term_id'] = adata.obs['Primary.Cardiac.Disease' if modality == 'atac' else 'Primary Cardiac Disease'].map(DISEASE_MAP)
    adata.obs['assay_ontology_term_id'] = adata.obs['sequencing_method'].map(ATAC_ASSAY_MAP if modality == 'atac' else RNA_ASSAY_MAP)
    adata.obs['self_reported_ethnicity_ontology_term_id'] = adata.obs['Ethnicity'].astype(str).map(ETHNICITY_MAP).astype('string').fillna('unknown').astype('category')    
    for col in adata.obs.columns[adata.obs.columns.str.endswith('_ontology_term_id')]:
        if adata.obs[col].isna().any():
            print(f"{title}: {col} has {adata.obs[col].isna().sum()} unmapped/NA values")
    adata.obs['DOI'] = adata.obs['Study'].map(DOI_MAP)
    print(f"{title}: unmapped studies:", adata.obs.loc[adata.obs['DOI'].isna(), 'Study'].unique())
    adata.uns['default_embedding'] = "X_umap"
    cell_col = 'vsmc_cell_type' if adata is vsmc_adata else 'vascular_cell_type'
    adata.obs['is_primary_data'] = False
    adata.uns['title'] = title
    adata.var['feature_is_filtered'] = False
    adata.uns['batch_condition'] = ['Study', 'donor_id', 'sample_type']

keep_obs = ['donor_id', 'Age', 'sex_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'BMI', 'disease_ontology_term_id', 'LVEF', 'MI_zone', 'development_stage_ontology_term_id', 'Study', 'DOI', 'tissue_ontology_term_id', 'tissue_type', 'sample_type', 'suspension_type', 'assay_ontology_term_id', 'multiome_GEX_linked', 'is_primary_data', 'cell_type_ontology_term_id', 'S_score', 'G2M_score']
keep_uns = ['title', 'organism_ontology_term_id', 'default_embedding', 'batch_condition']

for adata, title, modality in adatas:
    adata.obs.rename(columns={'linked_RNA': 'multiome_GEX_linked'}, inplace=True)
    cell_col = 'vsmc_cell_type' if adata is vsmc_adata else 'vascular_cell_type'
    adata.obs = adata.obs[[c for c in keep_obs + [cell_col] if c in adata.obs.columns]].copy()
    keep_obsm = {'X_umap', 'X_scVI' if adata is vsmc_adata else 'X_peakVI' if modality == 'atac' else 'X_scANVI'}
    for key in list(adata.obsm):
        if key not in keep_obsm: del adata.obsm[key]
    keep_uns_dataset = keep_uns + [f'{cell_col}_colors']
    for key in list(adata.uns):
        if key not in keep_uns_dataset: del adata.uns[key]
    adata.obsp.clear()
    adata.var = adata.var[[]].copy()
    adata.var['feature_is_filtered'] = False

# Change Colours to Match Manuscript
col = 'vascular_cell_type'
vascular_adata.obs[col] = vascular_adata.obs[col].astype('category')
rna_cats = list(vascular_adata.obs[col].cat.categories)
cmap = plt.get_cmap('tab20', len(rna_cats))
rna_color_map = {cat: mpl.colors.to_hex(cmap(i)) for i, cat in enumerate(rna_cats)}
vascular_adata.uns[f'{col}_colors'] = [rna_color_map[c] for c in rna_cats]
atac_color_map = {'Pericyte': 'Pericyte 1', 'Arterial Endothelial Cell': 'Arterial Endothelial Cell 1', 'Capillary Endothelial Cell': 'Capillary Endothelial Cell 2'}
atac_adata.obs[col] = atac_adata.obs[col].astype('category')
atac_adata.uns[f'{col}_colors'] = [rna_color_map[atac_color_map.get(c, c)] for c in atac_adata.obs[col].cat.categories]
vsmc_adata.obs['vsmc_cell_type'] = vsmc_adata.obs['vsmc_cell_type'].astype('category')
vsmc_adata.uns['vsmc_cell_type_colors'] = list(map(mpl.colors.to_hex, plt.get_cmap('Dark2').colors[:4]))

# Write Files to .h5ad
atac_adata.write_h5ad("/COMMA/results/data_submission/final/atac_CELLxGENE.h5ad")
vascular_adata.write_h5ad("/COMMA/results/data_submission/final/vascular_CELLxGENE.h5ad")
vsmc_adata.write_h5ad("/COMMA/results/data_submission/final/vsmc_CELLxGENE.h5ad")

# Filter Fragments File to Contain Only Relevant Barcodes
base = Path("/COMMA/data/meta_analysis/processed")
outdir = Path("/COMMA/results/data_submission/fragments")
outdir.mkdir(parents=True, exist_ok=True)

study_dirs = {
    "Amrute 2024": "multiome_amrute",
    "Kanemaru 2023": "multiome_kanemaru"}

meta = atac_adata.obs[["Study"]].copy()
meta["cell_id"] = atac_adata.obs_names
meta[["Sample", "barcode"]] = meta["cell_id"].str.rsplit("#", n=1, expand=True).to_numpy()

for (study, sample), cells in meta.groupby(["Study", "Sample"]):
    directory = base / study_dirs[study] / sample
    paths = [
        directory / "outs/atac_fragments.tsv.gz",
        directory / "outs/fragments.tsv.gz",
        directory / "fragments.tsv.gz"
    ]
    path = next((p for p in paths if p.exists()), None)
    if path is None:
        print(f"Missing: {study} {sample}")
        continue
    barcodes = dict(zip(cells["barcode"], cells["cell_id"]))
    output = outdir / f"{sample}.fragments.tsv"
    with gzip.open(path, "rt") as fin, open(output, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4 and fields[3] in barcodes:
                fields[3] = barcodes[fields[3]]
                _ = fout.write("\t".join(fields) + "\n")
    print(f"Saved: {output}")

# Concatenate Fragment Files
fragments = "/COMMA/results/data_submission/final/COMMA.fragments.tsv"
with open(fragments, "wb") as fout:
    for path in sorted(outdir.glob("*.fragments.tsv")):
        with open(path, "rb") as fin:
            shutil.copyfileobj(fin, fout)

# Check Fragments File is Complete Match
fragments = "/COMMA/results/data_submission/final/COMMA.fragments.tsv"
missing = set(atac_adata.obs_names.astype(str))
with open(fragments) as f:
    for line in f: missing.discard(line.split('\t')[3])

print(f"Matched: {atac_adata.n_obs - len(missing):,}/{atac_adata.n_obs:,} ({100 * (1 - len(missing)/atac_adata.n_obs):.2f}%)")
print("Unmatched:", list(missing)[:20])
extra = set()
obs = set(atac_adata.obs_names.astype(str))

with open(fragments) as f:
    for line in f:
        barcode = line.split('\t')[3]
        if barcode not in obs: extra.add(barcode)

print(f"Excess fragment barcodes: {len(extra):,}")
print("Examples:", list(extra)[:20])