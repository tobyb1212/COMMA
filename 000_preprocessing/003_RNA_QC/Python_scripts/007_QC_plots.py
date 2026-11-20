import pandas as pd
import anndata as ad
import scanpy as sc
import seaborn as sns
import os
import argparse
import random
import matplotlib.pyplot as plt

random.seed(1234)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process some directories.')
    parser.add_argument('-f', '--input_files', type=str, required=True, help='Comma-separated paths to input .h5ad Files')
    return parser.parse_args()

def save_plot(outdir, filename, dpi=300):
    os.makedirs(outdir, exist_ok=True)
    png_path = os.path.join(outdir, f"{filename}.png")
    pdf_path = os.path.join(outdir, f"{filename}.pdf")
    
    # Wrapper function to save files as .pdf and .png
    plt.savefig(png_path, dpi=dpi, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"Saved: {png_path} and {pdf_path}")

def main():
    args = parse_arguments()
    input_files = args.input_files.split(',')
    os.chdir("/COMMA/results/scrna/QC/")

    # Merge all study anndata objects together
    adata_list = [sc.read_h5ad(f) for f in input_files]
    combined_adata = adata_list[0].concatenate(*adata_list[1:], join='outer', index_unique=None)
    combined_adata.write("/COMMA/results/scrna/count_matrix/pre_filter_combined_count_matrix.h5ad")
    combined_adata = combined_adata[combined_adata.obs['Study']!= "Leblanc 2025"] # Generate plots without Leblanc as subsequently excluded

    # Plot metrics so they can be compared across studies here
    # Supplemental Figure 1E
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5), sharey=True, gridspec_kw={'width_ratios': [11, 3]})
    sc.pl.violin(combined_adata[combined_adata.obs['sample_type'] == "Nuclei"],["pct_counts_mt"], groupby='Study', stripplot=False, show=False, rotation=90, ax=axes[0])
    axes[0].set_ylim(0, 40) 
    axes[0].axhline(y=5, color='red', linestyle='dotted', linewidth=1)
    axes[0].set_title("Nuclei")
    
    sc.pl.violin(combined_adata[combined_adata.obs['sample_type'].isin(["CD45+ cells", "CD45- cells", "Cells"])],["pct_counts_mt"], groupby='Study', stripplot=False, show=False, rotation=90, ax=axes[1])
    axes[1].set_ylim(0, 40) 
    axes[1].axhline(y=20, color='red', linestyle='dotted', linewidth=1)
    axes[1].set_title("Cells")
    plt.tight_layout()
    save_plot(outdir=f'cell_metrics/', filename="mitochondrial_percentage")
    plt.close(fig)
    # Supplemental Figure 1E
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5), sharey=True,  gridspec_kw={'width_ratios': [11, 3]})
    sc.pl.violin(combined_adata[combined_adata.obs['sample_type'] == "Nuclei"],["pct_counts_ribo"], groupby='Study', stripplot=False, show=False, rotation=90, ax=axes[0])
    axes[0].set_ylim(0, 40) 
    axes[0].axhline(y=5, color='red', linestyle='dotted', linewidth=1)
    axes[0].set_title("Nuclei")
    
    sc.pl.violin(combined_adata[combined_adata.obs['sample_type'].isin(["CD45+ cells", "CD45- cells", "Cells"])],["pct_counts_ribo"], groupby='Study', stripplot=False, show=False, rotation=90, ax=axes[1])
    axes[1].set_ylim(0, 40) 
    axes[1].axhline(y=20, color='red', linestyle='dotted', linewidth=1)
    axes[1].set_title("Cells")
    plt.tight_layout()
    save_plot(outdir=f'cell_metrics/', filename="ribosomal_percentage")
    plt.close(fig)

    plt.figure(figsize=(7, 3.2))
    ax = sc.pl.violin(combined_adata,["pct_counts_hb"], groupby='Study', stripplot=False, show=False, rotation=90)
    ax.axhline(y=3, color='red', linestyle='dotted', linewidth=1)
    plt.ylim(0, 20) 
    save_plot(outdir=f'cell_metrics/', filename="haemoglobin_percentage")
    plt.close()

    # Supplemental Figure 1D
    plt.figure(figsize=(7, 3.2))
    ax = sc.pl.violin(combined_adata,["scrublet_score"], groupby='Study', stripplot=False, show=False, rotation=90)
    ax.axhline(y=0.3, color='red', linestyle='dotted', linewidth=1)
    plt.ylim(0, 0.8) 
    save_plot(outdir=f'cell_metrics/', filename="doublet_violin")
    plt.close()

    plt.figure(figsize=(7, 3.2))
    ax = sc.pl.violin(combined_adata,["total_counts"], stripplot=False, groupby='Study',show=False, rotation=90)
    ax.axhline(y=500, color='red', linestyle='dotted', linewidth=1)
    plt.ylim(0, 100000) 
    save_plot(outdir=f'cell_metrics/', filename="total_counts")
    plt.close()
    
    # Supplemental Figure 1C
    ax = sc.pl.violin(combined_adata,["stress_score"], groupby='Study', stripplot=False, show=False, rotation=90)
    plt.ylim(0, 20)
    save_plot(outdir=f'cell_metrics/', filename="stress_violin")
    plt.close()

    # Subset adata to only include unique fastq_names
    unique_mask = ~combined_adata.obs.duplicated(subset="fastq_name", keep="first")
    adata_unique_fastq = combined_adata[unique_mask].copy()

    # Plot violin plots on a per-study basis for key alignment metrics
    # Adding lines to indicate thresholds
    plt.figure(figsize=(7, 3.2))
    sc.pl.violin(adata_unique_fastq, keys='Number of Reads', groupby='Study', rotation=90)
    save_plot(outdir=f'alignment_metrics/', filename="number_of_reads")
    plt.close()

    plt.figure(figsize=(7, 3.2))
    sc.pl.violin(adata_unique_fastq, keys='Reads With Valid Barcodes', groupby='Study', rotation=90)
    save_plot(outdir=f'alignment_metrics/', filename="reads_valid_bc")
    plt.close()

    # Supplemental Figure 1A
    plt.figure(figsize=(7, 3.2))
    ax = sc.pl.violin(adata_unique_fastq, keys='Sequencing Saturation', groupby='Study', rotation=90, show=False)
    ax.axhline(y=0.2, color='red', linestyle='dotted', linewidth=1)
    save_plot(outdir=f'alignment_metrics/', filename="sequencing_saturation")
    plt.close()

    # Supplemental Figure 1B
    plt.figure(figsize=(7, 3.2))
    ax = sc.pl.violin(adata_unique_fastq, keys='Q30 Bases in RNA read', groupby='Study', rotation=90, show=False)
    ax.axhline(y=0.7, color='red', linestyle='dotted', linewidth=1)
    save_plot(outdir=f'alignment_metrics/', filename="q30_rna_read")
    plt.close()

    plt.figure(figsize=(7, 3.2))
    sc.pl.violin(adata_unique_fastq, keys='Reads Mapped to Genome: Unique', groupby='Study', rotation=90)
    save_plot(outdir=f'alignment_metrics/', filename="mapped_genome_unique")
    plt.close()

    plt.figure(figsize=(7, 3.2))
    ax = sc.pl.violin(adata_unique_fastq, keys='Reads Mapped to GeneFull_Ex50pAS: Unique GeneFull_Ex50pAS', groupby='Study', rotation=90, show=False)
    ax.axhline(y=0.3, color='red', linestyle='dotted', linewidth=1)
    save_plot(outdir=f'alignment_metrics/', filename="mapped_transcriptome_unique")
    plt.close()

    plt.figure(figsize=(7, 3.2))
    sc.pl.violin(adata_unique_fastq, keys='Estimated Number of Cells', groupby='Study', rotation=90)
    save_plot(outdir=f'alignment_metrics/', filename="estimated_cells")
    plt.close()

    plt.figure(figsize=(7, 3.2))
    sc.pl.violin(adata_unique_fastq, keys='Fraction of Unique Reads in Cells', groupby='Study', rotation=90)
    save_plot(outdir=f'alignment_metrics/', filename="fraction_reads_in_cells")
    plt.close()

    plt.figure(figsize=(7, 3.2))
    sc.pl.violin(adata_unique_fastq, keys='Mean Reads per Cell', groupby='Study', rotation=90)
    save_plot(outdir=f'alignment_metrics/', filename="reads_per_cell")
    plt.close()

    plt.figure(figsize=(7, 3.2))
    sc.pl.violin(adata_unique_fastq, keys='Median UMI per Cell', groupby='Study', rotation=90)
    save_plot(outdir=f'alignment_metrics/', filename="umi_per_cell")
    plt.close()

    plt.figure(figsize=(7, 3.2))
    sc.pl.violin(adata_unique_fastq, keys='Median GeneFull_Ex50pAS per Cell', groupby='Study', rotation=90)
    save_plot(outdir=f'alignment_metrics/', filename="median_gene_per_cell")
    plt.close()

    plt.figure(figsize=(7, 3.2))
    sc.pl.violin(adata_unique_fastq, keys='Total GeneFull_Ex50pAS Detected', groupby='Study', rotation=90)
    save_plot(outdir=f'alignment_metrics/', filename="total_genes")
    plt.close()

    print(f"Plotted all alignment metrics")
    print(f"Complete!")

if __name__ == "__main__":
    main()    