"""Render only panel 2 of the gene-embedding Euclidean distance change figure.

Input: All_Results workbook with pre/post_euclidean_distance_mean columns.
Output: histogram PDF/PNG and a compact JSON summary; no boxplot or scatter test.
Source: code_251216/20251211.py, distance-change histogram block.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from pipeline.common import write_json


def plot_distance_change(input_path, output_dir):
    table = pd.read_excel(input_path, sheet_name='All_Results')
    pre = table['pre_euclidean_distance_mean'].values
    post = table['post_euclidean_distance_mean'].values
    distance_change = pre - post
    if not len(distance_change) or not np.isfinite(distance_change).all():
        raise ValueError('Distance table must contain finite pre/post distances')
    decreased_count = np.sum(distance_change > 0)
    total_count = len(distance_change)
    decreased_ratio = decreased_count / total_count * 100
    fig, ax = plt.subplots(figsize=(3.5, 3.5))
    ax.hist(distance_change, bins=100, alpha=0.9, edgecolor='black', linewidth=1.0, color='#c15956')
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1.5, label='No change')
    ax.set_xlabel('Distance Change (Pre - Post)')
    ax.set_ylabel('Number of Genes')
    ax.set_title('Distribution of Distance Changes', fontsize=10)
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.get_major_formatter().set_scientific(False)
    ax.text(0.25, 0.95, f'Decreased: {decreased_count}/{total_count} ({decreased_ratio:.1f}%)\nMean: {np.mean(distance_change):.6f}',
            transform=ax.transAxes, va='top',
            bbox=dict(boxstyle='round', edgecolor='black', facecolor='white', alpha=0.9, linewidth=1), fontsize=8)
    fig.tight_layout()
    for extension in ['pdf', 'png']:
        fig.savefig(output_dir/f'euclidean_distance_change_panel2.{extension}', dpi=300, bbox_inches='tight')
    plt.close(fig)
    stats = {'n_genes':int(total_count), 'n_decreased':int(decreased_count),
             'decreased_percent':float(decreased_ratio), 'mean_change':float(np.mean(distance_change))}
    write_json(output_dir/'euclidean_distance_change_panel2.json', stats)
    return stats
