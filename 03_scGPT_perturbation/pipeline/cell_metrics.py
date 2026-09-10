"""Paired R/R-prime/NR metrics. Inputs: per-gene h5ad and gene output directories.
Outputs: all/filtered metrics and empirical comparisons. Use 05_kNN_enrichment.py.
The no-gate, NR gene-zero branch is linked by the selected-gene plotting script.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors



@dataclass
class PreparedGene:
    gene: str
    h5ad_path: Path
    status: str
    exclude_reason: str
    n_responder_nonzero: int
    n_responder_pert_nonzero: int
    n_responder_pert_for_test_nonzero: int
    n_nonresponder: int
    n_nonresponder_gene_zero: int
    n_paired_cells_tested: int
    emb_r: np.ndarray | None
    emb_rp: np.ndarray | None
    emb_nr: np.ndarray | None


def discover_gene_dirs(base_dir: Path) -> list[Path]:
    return sorted(p for p in base_dir.iterdir() if p.is_dir() and not p.name.startswith("__"))


def empirical_p_greater(observed: float, others: np.ndarray) -> float:
    return float((1 + np.sum(others >= observed)) / (len(others) + 1))

def rank_desc(observed: float, all_values: np.ndarray) -> int:
    return int(np.sum(all_values > observed) + 1)

def percentile_desc(observed: float, all_values: np.ndarray) -> float:
    return float(100.0 * np.sum(all_values <= observed) / len(all_values))

def cell_pair_key(cell_id: str) -> str:
    if "_" not in cell_id:
        return cell_id
    return cell_id.rsplit("_", 1)[0]

def find_gene_in_var(adata: ad.AnnData, gene: str) -> str | None:
    if gene in adata.var_names:
        return gene
    cand = [v for v in adata.var_names if str(v).lower() == gene.lower()]
    if len(cand) == 1:
        return cand[0]
    return None

def get_expr_vector(adata: ad.AnnData, gene: str) -> np.ndarray | None:
    g = find_gene_in_var(adata, gene)
    if g is None:
        return None
    x = adata[:, g].X
    if hasattr(x, "toarray"):
        x = x.toarray().ravel()
    else:
        x = np.asarray(x).ravel()
    return x

def prepare_gene(h5ad_path: Path, gene: str) -> PreparedGene:
    if not h5ad_path.exists():
        return PreparedGene(gene, h5ad_path, "missing", "h5ad_missing", 0, 0, 0, 0, 0, 0, None, None, None)

    try:
        adata = ad.read_h5ad(h5ad_path)
    except Exception as exc:
        return PreparedGene(gene, h5ad_path, "error", f"read_h5ad_error: {exc}", 0, 0, 0, 0, 0, 0, None, None, None)

    if "group" not in adata.obs.columns:
        return PreparedGene(gene, h5ad_path, "error", "missing_group", 0, 0, 0, 0, 0, 0, None, None, None)
    if "X_scGPT_cell" not in adata.obsm:
        return PreparedGene(gene, h5ad_path, "error", "missing_X_scGPT_cell", 0, 0, 0, 0, 0, 0, None, None, None)

    groups = adata.obs["group"].astype(str).to_numpy()
    cell_ids = adata.obs_names.astype(str).to_numpy()
    emb = adata.obsm["X_scGPT_cell"]
    if hasattr(emb, "toarray"):
        emb = emb.toarray()
    emb = np.asarray(emb)

    gene_expr = get_expr_vector(adata, gene)
    if gene_expr is None:
        return PreparedGene(gene, h5ad_path, "error", "perturbation_gene_not_found_in_var", 0, 0, 0, 0, 0, 0, None, None, None)

    nonzero_col = f"{gene}_nonzero"
    if nonzero_col in adata.obs.columns:
        nonzero_bool = adata.obs[nonzero_col].fillna(False).astype(bool).to_numpy()
    else:
        nonzero_bool = gene_expr > 0

    m_nr = groups == "nonresponder"
    m_r = groups == "responder"
    m_rp = groups == "responder_pert"
    m_nr_gene_zero = m_nr & (gene_expr == 0)
    m_rp_eligible = m_rp & nonzero_bool

    n_responder_nonzero = int(np.sum(m_r & nonzero_bool))
    n_responder_pert_nonzero = int(np.sum(m_rp & nonzero_bool))
    n_responder_pert_for_test_nonzero = int(np.sum(m_rp_eligible))
    n_nonresponder = int(np.sum(m_nr))
    n_nonresponder_gene_zero = int(np.sum(m_nr_gene_zero))

    responder_indices = np.where(m_r)[0]
    responder_pert_indices = np.where(m_rp_eligible)[0]
    responder_map = {cell_pair_key(cell_ids[i]): int(i) for i in responder_indices}
    responder_pert_map = {cell_pair_key(cell_ids[i]): int(i) for i in responder_pert_indices}
    common_keys = sorted(set(responder_map) & set(responder_pert_map))
    n_pairs = len(common_keys)

    if n_nonresponder_gene_zero <= 0:
        return PreparedGene(
            gene,
            h5ad_path,
            "error",
            "no_nonresponder_gene_zero_cells",
            n_responder_nonzero,
            n_responder_pert_nonzero,
            n_responder_pert_for_test_nonzero,
            n_nonresponder,
            n_nonresponder_gene_zero,
            n_pairs,
            None,
            None,
            None,
        )
    if n_pairs == 0:
        return PreparedGene(
            gene,
            h5ad_path,
            "error",
            "no_paired_responder_responder_pert_cells_after_gate",
            n_responder_nonzero,
            n_responder_pert_nonzero,
            n_responder_pert_for_test_nonzero,
            n_nonresponder,
            n_nonresponder_gene_zero,
            n_pairs,
            None,
            None,
            None,
        )

    responder_sel = np.array([responder_map[k] for k in common_keys], dtype=int)
    responder_pert_sel = np.array([responder_pert_map[k] for k in common_keys], dtype=int)
    nr_sel = np.where(m_nr_gene_zero)[0]

    return PreparedGene(
        gene,
        h5ad_path,
        "ok",
        "",
        n_responder_nonzero,
        n_responder_pert_nonzero,
        n_responder_pert_for_test_nonzero,
        n_nonresponder,
        n_nonresponder_gene_zero,
        n_pairs,
        emb[responder_sel],
        emb[responder_pert_sel],
        emb[nr_sel],
    )




def _nr_fraction_for_ref_queries(idx: np.ndarray, n_r: int, k: int) -> np.ndarray:
    out = np.zeros(idx.shape[0], dtype=float)
    for i in range(idx.shape[0]):
        neighbors = idx[i]
        neighbors = neighbors[neighbors != i][:k]
        if neighbors.size == 0:
            out[i] = np.nan
            continue
        out[i] = np.mean(neighbors >= n_r)
    return out

def _nr_fraction_for_nonref_queries(idx: np.ndarray, n_r: int) -> np.ndarray:
    return np.mean(idx >= n_r, axis=1)

def compute_knn_label_enrichment(prep: PreparedGene, k: int) -> dict[str, float]:
    out = {}
    ref = np.vstack([prep.emb_r, prep.emb_nr])
    n_r = prep.emb_r.shape[0]
    n_ref = ref.shape[0]
    if n_ref <= 1:
        out["knn_nr_fraction_shift_mean_euclidean"] = np.nan
        return out

    for metric in ("euclidean",):
        n_for_r = min(k + 1, n_ref)
        n_for_rp = min(k, n_ref)
        nn = NearestNeighbors(n_neighbors=max(n_for_r, n_for_rp), metric=metric)
        nn.fit(ref)

        idx_r = nn.kneighbors(prep.emb_r, n_neighbors=n_for_r, return_distance=False)
        frac_r = _nr_fraction_for_ref_queries(idx_r, n_r=n_r, k=min(k, n_ref - 1))

        idx_rp = nn.kneighbors(prep.emb_rp, n_neighbors=n_for_rp, return_distance=False)
        frac_rp = _nr_fraction_for_nonref_queries(idx_rp, n_r=n_r)

        valid = ~np.isnan(frac_r)
        if np.sum(valid) == 0:
            out[f"knn_nr_fraction_shift_mean_{metric}"] = np.nan
            continue

        shift = frac_rp[valid] - frac_r[valid]                              
        out[f"knn_nr_fraction_shift_mean_{metric}"] = float(np.mean(shift))

    return out


def run(args: argparse.Namespace) -> None:
    if not args.base_dir.exists():
        raise FileNotFoundError(f"Base dir not found: {args.base_dir}")
    if not args.input_dir.exists():
        raise FileNotFoundError(f"Input dir not found: {args.input_dir}")

    rows = []
    primary_metric: str | None = None

    for gdir in discover_gene_dirs(args.base_dir):
        gene = gdir.name
        h5ad_path = args.input_dir / f"{gene}_scGPT.h5ad"
        prep = prepare_gene(h5ad_path, gene)

        row: dict[str, object] = {
            "gene": gene,
            "h5ad_path": str(h5ad_path),
            "status": prep.status,
            "exclude_reason": prep.exclude_reason,
            "n_responder_nonzero": prep.n_responder_nonzero,
            "n_responder_pert_nonzero": prep.n_responder_pert_nonzero,
            "n_responder_pert_for_test_nonzero": prep.n_responder_pert_for_test_nonzero,
            "n_nonresponder": prep.n_nonresponder,
            "n_nonresponder_gene_zero": prep.n_nonresponder_gene_zero,
            "n_paired_cells_tested": prep.n_paired_cells_tested,
        }

        if prep.status == "ok" and prep.emb_r is not None and prep.emb_rp is not None and prep.emb_nr is not None:
            method_metrics = compute_knn_label_enrichment(prep, k=args.k)
            row.update(method_metrics)
            primary_metric = "knn_nr_fraction_shift_mean_euclidean"
        rows.append(row)

    all_df = pd.DataFrame(rows).sort_values("gene").reset_index(drop=True)

    all_df["pass_filter_n_responder_pert_nonzero_ge"] = (
        pd.to_numeric(all_df["n_responder_pert_for_test_nonzero"], errors="coerce") >= args.min_responder_pert_nonzero
    )
    all_df["pass_filter_n_nonresponder_gene_zero_gt"] = (
        pd.to_numeric(all_df["n_nonresponder_gene_zero"], errors="coerce") > args.min_nonresponder_gene_zero
    )
    all_df["is_target_gene"] = all_df["gene"] == args.target_gene

    if primary_metric is None:
        primary_metric = ""
    filtered_df = all_df[
        all_df["pass_filter_n_responder_pert_nonzero_ge"]
        & all_df["pass_filter_n_nonresponder_gene_zero_gt"]
        & (all_df[primary_metric].notna() if primary_metric in all_df.columns else False)
    ].copy()

    suffix = "knn_label_enrichment"
    tag = f"_{args.output_tag}" if args.output_tag else ""
    out_all_csv = args.base_dir / f"{suffix}_all_metrics_no_gate_nr_genezero_paired{tag}.csv"
    out_all_xlsx = args.base_dir / f"{suffix}_all_metrics_no_gate_nr_genezero_paired{tag}.xlsx"
    out_filt_csv = args.base_dir / f"{suffix}_filtered_metrics_no_gate_nr_genezero_paired{tag}.csv"
    out_filt_xlsx = args.base_dir / f"{suffix}_filtered_metrics_no_gate_nr_genezero_paired{tag}.xlsx"

    all_df.to_csv(out_all_csv, index=False)
    all_df.to_excel(out_all_xlsx, index=False, sheet_name="all_metrics")
    filtered_df.to_csv(out_filt_csv, index=False)
    filtered_df.to_excel(out_filt_xlsx, index=False, sheet_name="filtered_metrics")

    print(f"Saved: {out_all_csv}")
    print(f"Saved: {out_all_xlsx}")
    print(f"Saved: {out_filt_csv}")
    print(f"Saved: {out_filt_xlsx}")
    print("-" * 80)
