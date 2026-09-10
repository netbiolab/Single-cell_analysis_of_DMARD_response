"""Plot archived significance metrics; input XLSX, output PDF/PNG and summary CSV.
Public entry point: 08_plot_results.py. Original plotting/statistics functions retained.
"""

from __future__ import annotations

import argparse

from dataclasses import dataclass

from pathlib import Path

import matplotlib.pyplot as plt

import numpy as np

import pandas as pd

TARGET_COLORS = {
    "SERPINB2": "#1b9e77",
    "THBD": "#d95f02",
    "SIGLEC1": "#7570b3",
}

@dataclass(frozen=True)
class PlotSpec:
    input_path: Path
    metric: str
    prefix: str
    title: str
    x_label: str

def load_metric_table(path: Path, metric: str, target_genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    df = pd.read_excel(path)
    if "gene" not in df.columns:
        raise ValueError(f"'gene' column not found in {path}")
    if metric not in df.columns:
        raise ValueError(f"Metric column '{metric}' not found in {path}")

    work = df.loc[df[metric].notna(), ["gene", metric]].copy()
    work["gene"] = work["gene"].astype(str)
    work = work.sort_values(metric, ascending=False).reset_index(drop=True)
    work["rank_desc"] = np.arange(1, len(work) + 1)

    target_df = work.loc[work["gene"].isin(target_genes)].copy()
    missing = [g for g in target_genes if g not in set(target_df["gene"])]
    if missing:
        raise ValueError(f"Target genes not found in {path}: {', '.join(missing)}")

    target_df["gene"] = pd.Categorical(target_df["gene"], categories=target_genes, ordered=True)
    target_df = target_df.sort_values("gene").reset_index(drop=True)
    return work, target_df

def empirical_p_greater(observed: float, others: np.ndarray) -> float:
    return float((1.0 + np.sum(others >= observed)) / (len(others) + 1.0))

def compute_summaries(work: pd.DataFrame, target_df: pd.DataFrame, metric: str) -> pd.DataFrame:
    values = work[metric].to_numpy(dtype=float)
    rows: list[dict[str, object]] = []
    for _, target in target_df.iterrows():
        observed = float(target[metric])
        others = work.loc[work["gene"] != target["gene"], metric].to_numpy(dtype=float)
        rows.append(
            {
                "gene": str(target["gene"]),
                "metric_value": observed,
                "rank_desc": int(target["rank_desc"]),
                "n_genes": int(len(work)),
                "n_background": int(len(others)),
                "percentile_desc": float(100.0 * np.sum(values <= observed) / len(values)),
                "empirical_p_greater": empirical_p_greater(observed, others),
            }
        )
    return pd.DataFrame(rows)

def summary_text(summary_df: pd.DataFrame) -> str:
    lines = []
    for _, row in summary_df.iterrows():
        lines.append(
            f"{row['gene']}: value={row['metric_value']:.4g}, "
            f"rank={int(row['rank_desc'])}/{int(row['n_genes'])}, "
            f"pct={row['percentile_desc']:.1f}, p={row['empirical_p_greater']:.4g}"
        )
    return "\n".join(lines)

def draw_targets_on_rank(ax: plt.Axes, target_df: pd.DataFrame, metric: str) -> None:
    for _, row in target_df.iterrows():
        gene = str(row["gene"])
        color = TARGET_COLORS.get(gene, "#c62828")
        ax.scatter(row["rank_desc"], row[metric], s=58, color=color, edgecolors="black", linewidths=0.5, zorder=3)
        ax.annotate(
            gene,
            (row["rank_desc"], row[metric]),
            xytext=(6, 4),
            textcoords="offset points",
            fontsize=8.5,
            color=color,
            weight="bold",
        )

def stylize(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

def save_rank_plot(work: pd.DataFrame, target_df: pd.DataFrame, metric: str, spec: PlotSpec, out_dir: Path, summary_df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    ax.scatter(work["rank_desc"], work[metric], s=18, color="#c9c9c9", alpha=0.85, linewidths=0)
    draw_targets_on_rank(ax, target_df, metric)
    ax.set_xlabel("Rank (descending)")
    ax.set_ylabel(spec.x_label)
    ax.set_title(f"{spec.title}\nRanked distribution", fontsize=11)
    ax.text(
        0.98,
        0.98,
        summary_text(summary_df),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.6,
        bbox={"facecolor": "white", "edgecolor": "#d0d0d0", "boxstyle": "round,pad=0.3"},
    )
    stylize(ax)
    fig.tight_layout()
    fig.savefig(out_dir / f"{spec.prefix}_selected_ranked_plot.pdf")
    fig.savefig(out_dir / f"{spec.prefix}_selected_ranked_plot.png", dpi=300)
    plt.close(fig)

def save_summary_table(rows: list[dict[str, object]], out_dir: Path) -> None:
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "selected_genes_significance_plot_summary.csv", index=False)
    df.to_excel(out_dir / "selected_genes_significance_plot_summary.xlsx", index=False)
