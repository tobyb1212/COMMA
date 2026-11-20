import pandas as pd
import anndata as ad
import scanpy as sc
import seaborn as sns
import os
import argparse
import random
import matplotlib.pyplot as plt

random.seed(1234)

def save_plot(outdir, filename, dpi=300):
    os.makedirs(outdir, exist_ok=True)
    png_path = os.path.join(outdir, f"{filename}.png")
    pdf_path = os.path.join(outdir, f"{filename}.pdf")
    
    # Wrapper function to save files as .pdf and .png
    plt.savefig(png_path, dpi=dpi, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"Saved: {png_path} and {pdf_path}")

os.chdir("/COMMA/results/scrna/downstream_analysis/epicardial_coronary/QC/")
adata = sc.read_h5ad("/COMMA/results/scrna/downstream_analysis/epicardial_coronary/count_matrix/rna_count_matrix_metadata_scores.h5ad")

# Plot metrics so they can be compared across samples here
# Supplementary Figure 7E
plt.figure(figsize=(7, 3.2))
ax = sc.pl.violin(adata,["pct_counts_mt"], groupby='fastq_name', stripplot=False, show=False, rotation=90)
ax.axhline(y=20, color='red', linestyle='dotted', linewidth=1)
plt.ylim(0, 40) 
save_plot(outdir=f'cell_metrics/', filename="mitochondrial_percentage")
plt.close()

# Supplementary Figure 7D
plt.figure(figsize=(7, 3.2))
ax = sc.pl.violin(adata,["pct_counts_ribo"], groupby='fastq_name', stripplot=False, show=False, rotation=90)
ax.axhline(y=20, color='red', linestyle='dotted', linewidth=1)
plt.ylim(0, 40) 
save_plot(outdir=f'cell_metrics/', filename="ribosomal_percentage")
plt.close()

plt.figure(figsize=(7, 3.2))
ax = sc.pl.violin(adata,["pct_counts_hb"], groupby='fastq_name', stripplot=False, show=False, rotation=90)
ax.axhline(y=3, color='red', linestyle='dotted', linewidth=1)
plt.ylim(0, 20) 
save_plot(outdir=f'cell_metrics/', filename="haemoglobin_percentage")
plt.close()

# Supplementary Figure 7C
plt.figure(figsize=(7, 3.2))
ax = sc.pl.violin(adata,["scrublet_score"], groupby='fastq_name', stripplot=False, show=False, rotation=90)
ax.axhline(y=0.3, color='red', linestyle='dotted', linewidth=1)
plt.ylim(0, 0.8) 
save_plot(outdir=f'cell_metrics/', filename="doublet_violin")
plt.close()

plt.figure(figsize=(7, 3.2))
ax = sc.pl.violin(adata,["total_counts"], stripplot=False, groupby='fastq_name',show=False, rotation=90)
ax.axhline(y=500, color='red', linestyle='dotted', linewidth=1)
plt.ylim(0, 100000) 
save_plot(outdir=f'cell_metrics/', filename="total_counts")
plt.close()

# Supplementary Figure 7B
plt.figure(figsize=(7, 3.2))
ax = sc.pl.violin(adata,["stress_score"], stripplot=False, groupby='fastq_name',show=False, rotation=90)
ax.axhline(y=500, color='red', linestyle='dotted', linewidth=1)
plt.ylim(0, 100000) 
save_plot(outdir=f'cell_metrics/', filename="total_counts")
plt.close()

# Subset adata to only include unique fastq_names
unique_mask = ~adata.obs.duplicated(subset="fastq_name", keep="first")
adata_unique_fastq = adata[unique_mask].copy()

# Plot violin plots on a per-SAMPLE basis for key alignment metrics
# Adding lines to indicate thresholds
plt.figure(figsize=(7, 3.2))
sc.pl.violin(adata_unique_fastq, keys='Number of Reads', groupby='fastq_name', rotation=90)
save_plot(outdir=f'alignment_metrics/', filename="number_of_reads")
plt.close()

plt.figure(figsize=(7, 3.2))
sc.pl.violin(adata_unique_fastq, keys='Reads With Valid Barcodes', groupby='fastq_name', rotation=90)
save_plot(outdir=f'alignment_metrics/', filename="reads_valid_bc")
plt.close()

# Supplementary Figure 7A
plt.figure(figsize=(7, 3.2))
ax = sc.pl.violin(adata_unique_fastq, keys='Sequencing Saturation', groupby='fastq_name', rotation=90, show=False)
plt.ylim(0, 1)
ax.axhline(y=0.2, color='red', linestyle='dotted', linewidth=1)
save_plot(outdir=f'alignment_metrics/', filename="sequencing_saturation")
plt.close()
plt.figure(figsize=(7, 3.2))
ax = sc.pl.violin(adata_unique_fastq, keys='Q30 Bases in RNA read', groupby='fastq_name', rotation=90, show=False)
plt.ylim(0, 1)
ax.axhline(y=0.7, color='red', linestyle='dotted', linewidth=1)
save_plot(outdir=f'alignment_metrics/', filename="q30_rna_read")
plt.close()

plt.figure(figsize=(7, 3.2))
sc.pl.violin(adata_unique_fastq, keys='Reads Mapped to Genome: Unique', groupby='fastq_name', rotation=90)
save_plot(outdir=f'alignment_metrics/', filename="mapped_genome_unique")
plt.close()

plt.figure(figsize=(7, 3.2))
ax = sc.pl.violin(adata_unique_fastq, keys='Reads Mapped to GeneFull_Ex50pAS: Unique GeneFull_Ex50pAS', groupby='fastq_name', rotation=90, show=False)
save_plot(outdir=f'alignment_metrics/', filename="mapped_transcriptome_unique")
plt.close()

plt.figure(figsize=(7, 3.2))
sc.pl.violin(adata_unique_fastq, keys='Estimated Number of Cells', groupby='fastq_name', rotation=90)
save_plot(outdir=f'alignment_metrics/', filename="estimated_cells")
plt.close()

plt.figure(figsize=(7, 3.2))
sc.pl.violin(adata_unique_fastq, keys='Fraction of Unique Reads in Cells', groupby='fastq_name', rotation=90)
save_plot(outdir=f'alignment_metrics/', filename="fraction_reads_in_cells")
plt.close()

plt.figure(figsize=(7, 3.2))
sc.pl.violin(adata_unique_fastq, keys='Mean Reads per Cell', groupby='fastq_name', rotation=90)
save_plot(outdir=f'alignment_metrics/', filename="reads_per_cell")
plt.close()

plt.figure(figsize=(7, 3.2))
sc.pl.violin(adata_unique_fastq, keys='Median UMI per Cell', groupby='fastq_name', rotation=90)
save_plot(outdir=f'alignment_metrics/', filename="umi_per_cell")
plt.close()

plt.figure(figsize=(7, 3.2))
sc.pl.violin(adata_unique_fastq, keys='Median GeneFull_Ex50pAS per Cell', groupby='fastq_name', rotation=90)
save_plot(outdir=f'alignment_metrics/', filename="median_gene_per_cell")
plt.close()

plt.figure(figsize=(7, 3.2))
sc.pl.violin(adata_unique_fastq, keys='Total GeneFull_Ex50pAS Detected', groupby='fastq_name', rotation=90)
save_plot(outdir=f'alignment_metrics/', filename="total_genes")
plt.close()

    # Plot violin plots on a per-study basis for key cellbender metrics
plt.figure(figsize=(7, 3.2))
sc.pl.violin(adata_unique_fastq, keys='cellbender_fraction_reads_removed', groupby='fastq_name', rotation=90)
save_plot(outdir=f'cellbender_metrics/', filename="cellbender_fraction_reads_removed")
plt.close()

plt.figure(figsize=(7, 3.2))
sc.pl.violin(adata_unique_fastq, keys='cellbender_cells', groupby='fastq_name', rotation=90)
save_plot(outdir=f'cellbender_metrics/', filename="cellbender_total_cells")
plt.close()

plt.figure(figsize=(7, 3.2))
sc.pl.violin(adata_unique_fastq, keys='cellbender_fraction_non_empty_cells', groupby='fastq_name', rotation=90)
save_plot(outdir=f'cellbender_metrics/', filename="cellbender_fraction_nonempty_cells")
plt.close()
    
print(f"Plotted all alignment metrics")
print(f"Complete!")
