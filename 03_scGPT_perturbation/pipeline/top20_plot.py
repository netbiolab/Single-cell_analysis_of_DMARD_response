"""
Custom barplot for top 20 SPC euclidean mean genes
for selected perturbations: SIGLEC1, SERPINB2, and THBD.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

BASE_DIR = Path("outputs/SRC")
INPUT_XLSX = BASE_DIR / "spc_euclidean_top20_selected_genes.xlsx"
OUTPUT_PDF = BASE_DIR / "spc_euclidean_top20_selected_genes_barplot_custom.pdf"
OUTPUT_PNG = BASE_DIR / "spc_euclidean_top20_selected_genes_barplot_custom.png"

TARGET_GENES = ["SIGLEC1", "SERPINB2", "THBD"]
PANEL_COLORS = {
    "SIGLEC1": "#7570b3",
    "SERPINB2": "#1b9e77",
    "THBD": "#d95f02",
}

def extract_top20(df: pd.DataFrame, perturb_gene: str) -> tuple[list[str], list[float], float]:
    row_df = df.loc[df["perturb_gene"].astype(str) == perturb_gene]
    if row_df.empty:
        raise ValueError(f"Perturbation gene not found: {perturb_gene}")

    row = row_df.iloc[0]
    genes = [str(row[f"top{i}_gene"]) for i in range(1, 21)]
    values = [float(row[f"top{i}_spc_euclidean_mean"]) for i in range(1, 21)]
    summary = float(row["top20_mean_spc_euclidean_mean"])
    return genes, values, summary

def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_linewidth(1.5)
    ax.spines["top"].set_color("black")
    ax.spines["bottom"].set_linewidth(1.5)
    ax.spines["bottom"].set_color("black")
    ax.spines["left"].set_linewidth(1.5)
    ax.spines["left"].set_color("black")
    ax.spines["right"].set_visible(False)
    ax.set_facecolor("white")
    ax.tick_params(axis="x", labelsize=9)
    ax.tick_params(axis="y", length=0, labelsize=9)

def main() -> None:
    print("=== SPC Top 20 Selected Genes Barplot ===")
    print(f"Reading: {INPUT_XLSX}")

    df = pd.read_excel(INPUT_XLSX)

    fig, axes = plt.subplots(1, 3, figsize=(14, 7), sharex=False)

    global_max = 0.0
    extracted: list[tuple[str, list[str], list[float], float]] = []
    for perturb_gene in TARGET_GENES:
        genes, values, summary = extract_top20(df, perturb_gene)
        extracted.append((perturb_gene, genes, values, summary))
        global_max = max(global_max, max(values))

    x_limit = global_max * 1.08

    for ax, (perturb_gene, genes, values, summary) in zip(axes, extracted):
        color = PANEL_COLORS.get(perturb_gene, "#c55a5a")
        ax.barh(
            range(len(genes)),
            values,
            color=color,
            edgecolor="black",
            linewidth=1.1,
            alpha=0.9,
        )
        ax.set_xlim(0, x_limit)
        ax.set_yticks(range(len(genes)))
        ax.set_yticklabels(genes)
        ax.invert_yaxis()
        ax.set_title(
            f"{perturb_gene}\nTop20 mean = {summary:.4f}",
            fontsize=11,
            fontweight="bold",
            pad=10,
        )
        style_axis(ax)

    fig.text(
        0.5,
        0.04,
        "SPC Euclidean Mean",
        ha="center",
        fontsize=12,
        fontweight="bold",
    )
    fig.suptitle(
        "Top 20 Genes by SPC Euclidean Mean for Selected Perturbations",
        fontsize=14,
        fontweight="bold",
        y=0.98,
    )

    fig.tight_layout(rect=[0.02, 0.06, 1, 0.95], w_pad=2.0)
    fig.savefig(OUTPUT_PDF, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUTPUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(f"Saved PDF: {OUTPUT_PDF}")
    print(f"Saved PNG: {OUTPUT_PNG}")

if __name__ == "__main__":
    main()
