import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
import os
import argparse
import random

random.seed(1234)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process some directories.')
    parser.add_argument('-f', '--input_files', type=str, required=True, help='Comma-separated paths to input directories')
    parser.add_argument('-o', '--output_dir', required=True, help='Path to the output directory')
    parser.add_argument('-s', '--summaries', required=True, help='Comma-separated paths to summary CSV files')
    return parser.parse_args()

def extract_sample_id(file_path):
    """Extract sample_id from file path by identifying the part before 'Solo.out'"""
    path_parts = file_path.split(os.sep)  
    for part in path_parts:
        if "Solo.out" in part:  
            return part.split("Solo.out")[0]
    raise ValueError(f"Could not find 'Solo.out' in the path: {file_path}")

def load_summary_data(summary_csv_path):
    """Load summary CSV where metrics are rows rather than columns."""
    try:
        summary_df = pd.read_csv(summary_csv_path, header=None, index_col=0)  # Read without headers, using first column as index
        summary_dict = summary_df[1].to_dict()  # Convert second column to dictionary {metric_name: value}
        return summary_dict
    except Exception as e:
        print(f"Error reading summary file {summary_csv_path}: {e}")
        return {}

def get_common_summary_metrics(summary_csvs):
    """Find common row index names (metrics) across all summary CSV files."""
    all_metrics = [set(pd.read_csv(csv, header=None, index_col=0).index) for csv in summary_csvs]
    common_metrics = set.intersection(*all_metrics)  # Find common row names
    return common_metrics

def main():
    args = parse_arguments()
    input_files = args.input_files.split(',')
    summary_csvs = args.summaries.split(',')
    output_dir = args.output_dir

    if len(input_files) != len(summary_csvs):
        raise ValueError("Mismatch between number of input files and summary CSVs!")

    # Find common metrics in all summary files
    common_metrics = get_common_summary_metrics(summary_csvs)
    print(f"Common metrics in summary files: {common_metrics}")

    adatas = {}

    for filename, summary_csv in zip(input_files, summary_csvs):
        sample_id = extract_sample_id(filename)
        print(f"Processing: {sample_id}")

        # Load cellbender counts into adata
        sample_adata = sc.read_h5ad(filename)

        # Load raw gzipped STARsolo counts and add as a layer to adata
        parent_dir = os.path.dirname(os.path.abspath(filename))
        input_folder = os.path.join(parent_dir, "../raw_gzip")
        sample_adata_star = sc.read_10x_mtx(input_folder)
        sample_adata.layers["starsolo"] = sample_adata_star[sample_adata.obs_names, :].X

        # Change var and obs names
        sample_adata.var_names_make_unique()
        sample_adata.obs.index = [f"{sample_id}_{bc}" for bc in sample_adata.obs.index]

        # Load summary metadata and filter to only contain summary metrics present in all files
        summary_data = load_summary_data(summary_csv)
        filtered_summary_data = {k: v for k, v in summary_data.items() if k in common_metrics}

        # Add metadata to AnnData object
        if filtered_summary_data:
            for key, value in filtered_summary_data.items():
                sample_adata.obs[key] = value  # Apply metadata to all cells from this sample

        # Add cellbender metrics to adata.obs
        cellbender_path = os.path.join(os.path.dirname(filename), 'cellbender_output_metrics.csv')
        cellbender_metrics = pd.read_csv(cellbender_path)
        fraction_reads_removed = cellbender_metrics.iloc[3, 1] 
        found_cells = cellbender_metrics.iloc[10, 1] 
        fraction_nonempty_cells = cellbender_metrics.iloc[14, 1]
        sample_adata.obs['cellbender_fraction_reads_removed'] = fraction_reads_removed
        sample_adata.obs['cellbender_cells'] = found_cells
        sample_adata.obs['cellbender_fraction_non_empty_cells'] = fraction_nonempty_cells

        adatas[sample_id] = sample_adata

    # Merge all outputs together
    adata = ad.concat(adatas, label="fastq_name")
    adata.obs_names_make_unique()
    print(adata.obs["fastq_name"].value_counts())

    # Compute basic QC metrics & filter file
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]", regex=True)
    
    sc.pp.filter_cells(adata, min_genes=200, inplace=True)
    sc.pp.filter_cells(adata, min_counts=500, inplace=True)
    sc.pp.filter_genes(adata, min_cells=50, inplace=True)

    # Save processed data to .h5ad file
    print(adata.shape)
    adata.write(os.path.join(output_dir, 'rna_count_matrix.h5ad'))
    print(f"Processed data saved to {output_dir}/rna_count_matrix.h5ad")

if __name__ == "__main__":
    main()
