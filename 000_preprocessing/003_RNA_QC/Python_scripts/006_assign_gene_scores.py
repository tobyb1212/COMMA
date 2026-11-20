import pandas as pd
import anndata
import scanpy as sc
import os
import random
import argparse

random.seed(1234)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process some directories.')
    parser.add_argument('-f', '--input_file', type=str, required=True, help='Individual study anndata object')
    return parser.parse_args()

def main():
    args = parse_arguments()
    input_file = args.input_file

    # Read anndata input file
    adata = sc.read_h5ad(input_file)

    # Score G2M & S Phase using Cell Cycle Genes from Seurat
    g2m_genes = pd.read_csv("/COMMA/results/scrna/metadata/g2m_genes.csv", header=0).iloc[:, 0].dropna().tolist()
    s_genes = pd.read_csv("/COMMA/results/scrna/metadata/s_genes.csv", header=0).iloc[:, 0].dropna().tolist()
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    
    # Score stress-related genes, gene set from O'Flanagan 2019, https://doi.org/10.1186/s13059-019-1830-0
    # .csv file from https://github.com/kieranrcampbell/scrnaseq-digestion-paper/blob/master/data/deliverables/coregene_df-FALSE-v3.csv
    stress_genes=pd.read_csv("/COMMA/results/scrna/metadata/coregene_df-FALSE-v3.csv", header=0).iloc[:, 6].dropna().tolist()
    sc.tl.score_genes(adata, gene_list=stress_genes, ctrl_size=len(stress_genes), score_name='stress_score', random_state = 1234)

    # Calculate Quality Control Metrics for Percentages of Mitochondrial, Ribosomal and Haemoglobin Genes
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

    # Write count matrix to .h5ad of the same name
    adata.write_h5ad(input_file)
    print(f"Processed file: {input_file}")

if __name__ == "__main__":
    main()