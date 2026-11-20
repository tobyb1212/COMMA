import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
import seaborn as sns
import os
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process scib outputs.')
    parser.add_argument('-f', '--input_csv', type=str, required=True, help='Comma-separated paths to input csv file outputs from scib')
    return parser.parse_args()

# Categorise all metrics as either Biological Conservation or Batch Effect Removal
BIO_METRICS = {'NMI', 'ARI', 'ASW (Labels)', 'Cell Cycle Conservation', 'HVG Conservation', 'Isolated Labels F1', 'Isolated Labels ASW', 'Trajectory Conservation', 'cLISI'}
TECH_METRICS = {'kBET', 'Graph Connectivity','ASW (Batch)', 'iLISI', 'PCR'}
def categorize_metric(metric_name):
    if metric_name in BIO_METRICS:
        return 'Biological Conservation'
    elif metric_name in TECH_METRICS:
        return 'Batch Correction'
    else:
        return 'Unknown'

# Rename Metrics to Correct Name
METRIC_RENAME_MAP = {
    'NMI_cluster/label': 'NMI',
    'ARI_cluster/label': 'ARI',
    'ASW_label': 'ASW (Labels)',
    'cell_cycle_conservation': 'Cell Cycle Conservation',
    'cLISI': 'cLISI',
    'kBET': 'kBET',
    'graph_conn': 'Graph Connectivity',
    'ASW_label/batch': 'ASW (Batch)',
    'iLISI': 'iLISI',
    'PCR_batch': 'PCR',}
def rename_metric(metric_name):
    return METRIC_RENAME_MAP.get(metric_name, metric_name)

# Wrap labels to look nice on graph
def wrap_label(label, max_len=14):
    words = label.split()
    lines = []
    current_line = ""
    for word in words:
        if len(current_line + ' ' + word) <= max_len if current_line else len(word) <= max_len:
            current_line = (current_line + ' ' + word).strip()
        else:
            lines.append(current_line)
            current_line = word
    if current_line:
        lines.append(current_line)
    return '\n'.join(lines)

# Presentation of this data is inspired by https://github.com/YosefLab/scib-metrics
def plot_results_table(merged_df: pd.DataFrame, summary_df: pd.DataFrame, save_path: str = None):
    import matplotlib.pyplot as plt
    from matplotlib.colorbar import ColorbarBase
    from matplotlib.colors import Normalize

    plot_df = merged_df.pivot(index='Method', columns='Metric', values='Score').astype(float)

    bio_order = [m for m in BIO_METRICS if m in plot_df.columns]
    tech_order = [m for m in TECH_METRICS if m in plot_df.columns]
    score_cols = ['Biological Conservation', 'Batch Correction', 'Total']
    plot_df = plot_df[bio_order + tech_order]

    methods = plot_df.index.tolist()
    metrics = plot_df.columns.tolist()
    total_width = len(metrics) + len(score_cols)

    # Colormaps
    cmap = plt.cm.PRGn
    norm = Normalize(vmin=0, vmax=1)
    cmap_bar = plt.cm.YlGnBu_r

    # Figure
    fig, ax = plt.subplots(figsize=(1.3 * total_width, 0.7 * len(methods) + 2))
    ax.set_xlim(-2, total_width - 0.5)
    ax.set_ylim(-1.5, len(methods) - 0.5)
    ax.invert_yaxis()
    ax.axis('off')

    # Y-axis (method) labels — centered, horizontal, multiline if needed
    for i, method in enumerate(methods):
        wrapped_label = '\n'.join(method.split())  # optional word wrap
        ax.text(-2.5, i, wrapped_label, va='center', ha='left', fontsize=10, weight='bold', rotation=0)

    # Metric score bubbles
    for i, method in enumerate(methods):
        for j, metric in enumerate(metrics):
            score = plot_df.loc[method, metric]
            radius = 0.25
            color = cmap(norm(score))
            ax.add_patch(plt.Circle((j, i), radius=radius, color=color))
            ax.text(j, i, f"{score:.2f}", ha='center', va='center', fontsize=8,
                    color='black' if 0.25 <= score <= 0.75 else 'white')

    # Bars for aggregate scores
    for k, col in enumerate(score_cols):
        col_idx = len(metrics) + k
        for i, method in enumerate(methods):
            score = summary_df.loc[summary_df['Method'] == method, col].values[0]
            bar_left = col_idx - 0.4  # shift left slightly to leave space
            ax.barh(i, score, left=bar_left, height=0.4, color=cmap_bar(score), align='center')
            ax.text(bar_left + 0.05, i, f"{score:.2f}", va='center', ha='left', fontsize=8)

    # Row dividers
    for i in range(len(methods)):
        ax.axhline(i - 0.5, color='gray', linestyle=':', linewidth=0.5)

    # Vertical dividers
    ax.axvline(len(bio_order) - 0.5, color='black', linestyle='-', linewidth=1)
    ax.axvline(len(metrics) - 0.5, color='black', linestyle='-', linewidth=1)

    # Section headers
    ax.text(len(bio_order) / 2 - 0.5, -1.5, 'Biological Conservation', ha='center', fontsize=13, weight='bold')
    ax.text(len(bio_order) + len(tech_order) / 2 - 0.5, -1.5, 'Batch Correction', ha='center', fontsize=13, weight='bold')
    ax.text(len(metrics) + len(score_cols) / 2 - 0.5, -1.5, 'Aggregate Score', ha='center', fontsize=13, weight='bold')

    # Top axis labels (horizontal)
    for j, label in enumerate(metrics + score_cols):
        wrapped = wrap_label(label)
        ax.text(j,-0.7, wrapped, rotation=0, ha='center', va='bottom', fontsize=8)

    ax.axhline(-0.5, color='black', linestyle='-', linewidth=1)

    # Colorbar (tiny, bottom left)
    cbar_ax = fig.add_axes([0.2, 0.02, 0.08, 0.02])  # [left, bottom, width, height]
    cb = ColorbarBase(cbar_ax, cmap=cmap, norm=norm, orientation='horizontal')
    cb.outline.set_linewidth(0.5)
    cbar_ax.tick_params(labelsize=6, length=2)
    cbar_ax.set_xticks([0, 0.5, 1])
    cbar_ax.set_xticklabels(["0", "0.5", "1"])

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()

def main():
    args = parse_arguments()
    input_files = args.input_csv.split(',')

    all_metrics = []

    # Read individual metrics file and extract the method
    for filename in input_files:
        df = pd.read_csv(filename, index_col=0)
        df = df.reset_index()
        df.columns = ['Metric', 'Score']
        df = df.dropna(subset=['Score'])
        method = '_'.join(os.path.basename(filename).split('_')[:-1])
        df['Method'] = method
        df['Metric'] = df['Metric'].apply(rename_metric)
        df = df[df['Metric'] != 'PCR']
        df['Category'] = df['Metric'].apply(categorize_metric)
        all_metrics.append(df)
        
    merged_df = pd.concat(all_metrics, ignore_index=True)

    summary = merged_df.groupby(['Method', 'Category'])['Score'].mean().unstack()
    summary['Total'] = summary.mean(axis=1)
    summary = summary.reset_index()

    # Supplemental Figure 3
    plot_results_table(merged_df, summary, save_path="scib/scib_results.png")
    plot_results_table(merged_df, summary, save_path="scib/scib_results.pdf")

if __name__ == "__main__":
    main()