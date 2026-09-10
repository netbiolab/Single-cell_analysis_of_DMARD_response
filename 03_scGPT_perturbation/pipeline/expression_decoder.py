"""Train expression predictor from scGPT embeddings and predict perturbed cells.

INPUT: paired AnnData with obs[group] and expression targets; CLS + gene embeddings.
OUTPUT: predicted R-prime expression H5AD/NPY/PT, training metrics and checkpoint.

This script trains on non-perturbed cells (group != responder_pert) and predicts
expression for perturbed cells (group == responder_pert) using CLS + gene embedding
features with a shared gene-wise MLP decoder.
"""
from __future__ import annotations
import argparse
import csv
import json
import logging
import random
import sys
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
import torch.nn as nn
from torch.cuda.amp import GradScaler, autocast
from torch.utils.data import DataLoader, Dataset
LOGGER = logging.getLogger('train_predict_expression')

def configure_logging() -> None:
    """Configure console logger."""
    logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s', handlers=[logging.StreamHandler(sys.stdout)])

def set_seed(seed: int) -> None:
    """Set random seeds for reproducibility."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

def _normalize_gene_name_list(obj: Any) -> Optional[List[str]]:
    """Normalize possible gene metadata into a list of strings."""
    if obj is None:
        return None
    if isinstance(obj, torch.Tensor):
        if obj.ndim != 1:
            return None
        try:
            obj = obj.cpu().numpy()
        except Exception:
            return None
    if isinstance(obj, np.ndarray):
        if obj.ndim != 1:
            return None
        return [str(x) for x in obj.tolist()]
    if isinstance(obj, pd.Index):
        return [str(x) for x in obj.tolist()]
    if isinstance(obj, (list, tuple)):
        return [str(x) for x in obj]
    return None

def load_embedding_tensor(embedding_path: str) -> Tuple[torch.Tensor, Optional[List[str]], Dict[str, Any]]:
    """Load embedding tensor from .pt and optionally collect gene metadata.

    Returns:
        emb_tensor: shape (n_cells, n_genes_plus_cls, emb_dim)
        emb_gene_names: optional gene name list for embedding gene axis (CLS excluded if needed)
        metadata: lightweight info dict for debugging
    """
    path = Path(embedding_path)
    if not path.exists():
        raise FileNotFoundError(f'Embedding file not found: {embedding_path}')
    raw = torch.from_numpy(np.load(path, mmap_mode='c', allow_pickle=False)) if path.suffix == '.npy' else torch.load(path, map_location='cpu')
    metadata: Dict[str, Any] = {'source_type': str(type(raw))}
    emb: Optional[torch.Tensor] = None
    gene_names: Optional[List[str]] = None
    if isinstance(raw, torch.Tensor):
        emb = raw
        metadata['tensor_key'] = None
    elif isinstance(raw, dict):
        keys = list(raw.keys())
        metadata['dict_keys'] = keys
        tensor_key_candidates = ['gene_emb', 'embeddings', 'embedding', 'emb', 'gene_embeddings', 'output', 'tensor', 'x']
        found_tensor_key: Optional[str] = None
        for key in tensor_key_candidates:
            val = raw.get(key)
            if isinstance(val, torch.Tensor) and val.ndim == 3:
                found_tensor_key = key
                emb = val
                break
        if emb is None:
            for key, val in raw.items():
                if isinstance(val, torch.Tensor) and val.ndim == 3:
                    found_tensor_key = key
                    emb = val
                    break
        metadata['tensor_key'] = found_tensor_key
        if emb is None:
            raise ValueError(f'Failed to extract 3D embedding tensor from dict .pt. Available keys: {keys}')
        gene_key_candidates = ['gene_names', 'var_names', 'genes', 'gene_symbols', 'feature_names', 'tokens', 'vocab']
        for key in gene_key_candidates:
            names = _normalize_gene_name_list(raw.get(key))
            if names is not None:
                gene_names = names
                metadata['gene_name_key'] = key
                break
    else:
        raise TypeError(f'Unsupported type loaded from {embedding_path}: {type(raw)}. Expected torch.Tensor or dict.')
    assert emb is not None
    if emb.ndim != 3:
        raise ValueError(f'Embedding tensor must be 3D (cell, gene+cls, emb_dim), got shape={tuple(emb.shape)}')
    emb = emb.float().contiguous()
    LOGGER.info('[load] embedding shape=%s, dtype=%s, source=%s', tuple(emb.shape), emb.dtype, metadata.get('source_type'))
    if gene_names is not None:
        LOGGER.info('[load] detected embedding gene metadata: n=%d', len(gene_names))
    return (emb, gene_names, metadata)

def _strip_possible_cls(names: List[str], expected_gene_count: int) -> Optional[List[str]]:
    """Return gene names aligned to gene axis (CLS excluded when necessary)."""
    if len(names) == expected_gene_count:
        return names
    if len(names) == expected_gene_count + 1:
        return names[1:]
    return None

def validate_alignment(adata: ad.AnnData, emb_tensor: torch.Tensor, emb_gene_names: Optional[List[str]]) -> Tuple[torch.Tensor, np.ndarray, str]:
    """Validate cell/gene alignment and return aligned embedding tensor.

    Strategy:
      1) strict cell match is mandatory
      2) if gene counts match, use direct adata.var order assumption
      3) if mismatch, try exact 1:1 mapping with embedding gene metadata
      4) otherwise fail with explicit error

    Returns:
        aligned_emb: embedding tensor aligned to adata.var order
        gene_indices_in_emb: indices used on embedding gene axis (excluding CLS)
        alignment_note: human-readable alignment assumption
    """
    n_cells_emb, n_gene_slots_emb, _ = emb_tensor.shape
    n_genes_emb = n_gene_slots_emb - 1
    LOGGER.info('[validate] adata shape=%s', adata.shape)
    LOGGER.info('[validate] embedding shape=%s (gene_excl_cls=%d)', tuple(emb_tensor.shape), n_genes_emb)
    if n_cells_emb != adata.n_obs:
        raise ValueError(f'Cell count mismatch. adata.n_obs={adata.n_obs}, embedding_cells={n_cells_emb}. Embedding must be extracted from the same adata ordering.')
    if n_genes_emb == adata.n_vars:
        note = 'Embedding gene axis (excluding CLS) is assumed to match adata.var_names order exactly.'
        LOGGER.info('[validate] gene count match: %d', n_genes_emb)
        LOGGER.info('[validate] alignment assumption: %s', note)
        return (emb_tensor, np.arange(n_genes_emb, dtype=np.int64), note)
    LOGGER.error('[validate] gene count mismatch: adata.n_vars=%d vs emb_genes=%d', adata.n_vars, n_genes_emb)
    if emb_gene_names is None:
        raise ValueError('Gene count mismatch and embedding gene metadata not found in .pt. Cannot infer mapping safely. Aborting.')
    stripped = _strip_possible_cls(emb_gene_names, n_genes_emb)
    if stripped is None:
        raise ValueError(f'Embedding gene metadata length is inconsistent with embedding tensor. len(gene_names)={len(emb_gene_names)}, expected {n_genes_emb} or {n_genes_emb + 1}.')
    emb_names = np.asarray(stripped, dtype=object)
    adata_names = np.asarray(adata.var_names.tolist(), dtype=object)
    uniq_emb, counts = np.unique(emb_names, return_counts=True)
    if (counts > 1).any():
        dup = uniq_emb[counts > 1][:10].tolist()
        raise ValueError(f'Embedding gene metadata contains duplicates; exact mapping impossible. Example duplicates={dup}')
    idx_map = {g: i for i, g in enumerate(emb_names.tolist())}
    missing = [g for g in adata_names.tolist() if g not in idx_map]
    if missing:
        raise ValueError(f'Cannot build complete gene mapping from embedding metadata to adata.var_names. Missing genes={len(missing)} (examples={missing[:10]}).')
    mapped = np.array([idx_map[g] for g in adata_names.tolist()], dtype=np.int64)
    if np.unique(mapped).size != adata.n_vars:
        raise ValueError('Mapped embedding gene indices are not unique. Only strict 1:1 mapping is allowed.')
    LOGGER.warning('[validate] applying explicit gene-name mapping from embedding metadata.')
    gene_part = emb_tensor[:, 1:, :]
    aligned_gene_part = gene_part[:, mapped, :]
    aligned_emb = torch.cat([emb_tensor[:, :1, :], aligned_gene_part], dim=1).contiguous()
    note = 'Embedding gene axis was reordered using exact gene-name mapping to adata.var_names. Mapping required complete 1:1 coverage and uniqueness.'
    LOGGER.info('[validate] alignment assumption: %s', note)
    return (aligned_emb, mapped, note)

def _select_target_source(adata: ad.AnnData, target_mode: str) -> Tuple[Any, str]:
    """Select target matrix source from AnnData according to target_mode."""
    if 'data' in adata.layers:
        return (adata.layers['data'], "layers['data']")
    return (adata.X, 'X')

def _row_to_numpy(row: Any, dtype: np.dtype) -> np.ndarray:
    """Convert one target row to dense 1D numpy array with desired dtype."""
    if sp.issparse(row):
        arr = row.toarray()
    else:
        arr = np.asarray(row)
    arr = np.asarray(arr).reshape(-1)
    return arr.astype(dtype, copy=False)

def estimate_dense_bytes(n_rows: int, n_cols: int, dtype: np.dtype) -> int:
    """Estimate memory bytes for dense matrix."""
    return int(n_rows) * int(n_cols) * np.dtype(dtype).itemsize

def split_indices(indices: np.ndarray, labels: np.ndarray, train_ratio: float, seed: int, stratified: bool) -> Tuple[np.ndarray, np.ndarray]:
    """Split indices into train/val with optional stratification."""
    rng = np.random.RandomState(seed)
    indices = np.asarray(indices, dtype=np.int64)
    labels = np.asarray(labels)
    val_ratio = 1.0 - train_ratio
    if not stratified:
        perm = indices.copy()
        rng.shuffle(perm)
        n_val = max(1, int(round(len(perm) * val_ratio)))
        n_val = min(n_val, len(perm) - 1)
        val_idx = np.sort(perm[:n_val])
        train_idx = np.sort(perm[n_val:])
        return (train_idx, val_idx)
    train_parts: List[np.ndarray] = []
    val_parts: List[np.ndarray] = []
    for group in np.unique(labels):
        group_idx = indices[labels == group]
        local = group_idx.copy()
        rng.shuffle(local)
        n_val_group = max(1, int(round(len(local) * val_ratio)))
        if n_val_group >= len(local):
            n_val_group = len(local) - 1
        if n_val_group <= 0:
            raise ValueError(f'Stratified split not possible for group={group}.')
        val_parts.append(local[:n_val_group])
        train_parts.append(local[n_val_group:])
    train_idx = np.sort(np.concatenate(train_parts))
    val_idx = np.sort(np.concatenate(val_parts))
    return (train_idx, val_idx)

class ExpressionDataset(Dataset):
    """Dataset for embedding tensor and target matrix with on-the-fly row densification."""

    def __init__(self, embeddings: torch.Tensor, target_source: Any, cell_indices: np.ndarray, n_genes: int, target_mode: str) -> None:
        self.embeddings = embeddings
        self.target_source = target_source
        self.cell_indices = np.asarray(cell_indices, dtype=np.int64)
        self.n_genes = int(n_genes)
        self.target_mode = target_mode
        self.target_dtype = np.float32
        if self.embeddings.shape[0] != len(self.cell_indices):
            raise ValueError(f'Embedding subset size and cell_indices length mismatch: {self.embeddings.shape[0]} vs {len(self.cell_indices)}')
        if self.embeddings.shape[1] - 1 != self.n_genes:
            raise ValueError(f'Embedding gene count mismatch with target genes: emb_genes={self.embeddings.shape[1] - 1}, target_genes={self.n_genes}')

    def __len__(self) -> int:
        return self.embeddings.shape[0]

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        emb_row = self.embeddings[idx]
        cell_idx = int(self.cell_indices[idx])
        row = self.target_source[cell_idx]
        target_np = _row_to_numpy(row, self.target_dtype)
        if target_np.shape[0] != self.n_genes:
            raise ValueError(f'Target row gene count mismatch at cell={cell_idx}: {target_np.shape[0]} vs expected {self.n_genes}')
        target_t = torch.from_numpy(target_np.astype(np.float32, copy=False))
        return (emb_row, target_t)

class SharedGeneDecoder(nn.Module):
    """Shared MLP decoder applied identically across all genes."""

    def __init__(self, feature_dim: int, hidden_dim: int, output_dim: int, num_layers: int=2, dropout: float=0.1) -> None:
        super().__init__()
        if num_layers < 1:
            raise ValueError('num_layers must be >= 1')
        layers: List[nn.Module] = [nn.LayerNorm(feature_dim), nn.Linear(feature_dim, hidden_dim), nn.GELU(), nn.Dropout(dropout)]
        for _ in range(num_layers - 1):
            layers.extend([nn.Linear(hidden_dim, hidden_dim), nn.GELU(), nn.Dropout(dropout)])
        layers.append(nn.Linear(hidden_dim, output_dim))
        self.mlp = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Forward for x shape (B, G, F), output (B, G, O)."""
        bsz, n_genes, feat_dim = x.shape
        y = self.mlp(x.reshape(bsz * n_genes, feat_dim))
        return y.reshape(bsz, n_genes, -1)

class ExpressionPredictor(nn.Module):
    """CLS + gene embedding expression predictor with shared gene decoder."""

    def __init__(self, emb_dim: int, hidden_dim: int, target_mode: str, n_bins: Optional[int], num_layers: int, dropout: float) -> None:
        super().__init__()
        self.emb_dim = emb_dim
        self.target_mode = target_mode
        feature_dim = 3 * emb_dim
        output_dim = 1
        self.decoder = SharedGeneDecoder(feature_dim=feature_dim, hidden_dim=hidden_dim, output_dim=output_dim, num_layers=num_layers, dropout=dropout)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Forward.

        Args:
            x: (B, G+1, D)
        Returns:
            
            data: (B, G) regression output
        """
        if x.ndim != 3:
            raise ValueError(f'Input must be 3D, got shape={tuple(x.shape)}')
        _, g_plus_cls, d = x.shape
        if d != self.emb_dim:
            raise ValueError(f'Embedding dim mismatch: model={self.emb_dim}, input={d}')
        cls = x[:, 0, :]
        gene = x[:, 1:, :]
        n_genes = g_plus_cls - 1
        cls_expanded = cls.unsqueeze(1).expand(-1, n_genes, -1)
        feat = torch.cat([gene, cls_expanded, gene * cls_expanded], dim=-1)
        out = self.decoder(feat)
        return out.squeeze(-1)

def _compute_loss(pred: torch.Tensor, target: torch.Tensor, loss_fn: nn.Module, target_mode: str) -> torch.Tensor:
    """Compute loss with mode-specific tensor shape handling."""
    return loss_fn(pred, target)

def train_one_epoch(model: nn.Module, loader: DataLoader, optimizer: torch.optim.Optimizer, loss_fn: nn.Module, device: torch.device, target_mode: str, scaler: Optional[GradScaler], grad_clip: float) -> float:
    """Train one epoch and return average loss."""
    model.train()
    total_loss = 0.0
    steps = 0
    for emb, target in loader:
        emb = emb.to(device, non_blocking=True)
        target = target.to(device, non_blocking=True)
        optimizer.zero_grad(set_to_none=True)
        if scaler is not None:
            with autocast():
                pred = model(emb)
                loss = _compute_loss(pred, target, loss_fn, target_mode)
            scaler.scale(loss).backward()
            if grad_clip > 0:
                scaler.unscale_(optimizer)
                nn.utils.clip_grad_norm_(model.parameters(), max_norm=grad_clip)
            scaler.step(optimizer)
            scaler.update()
        else:
            pred = model(emb)
            loss = _compute_loss(pred, target, loss_fn, target_mode)
            loss.backward()
            if grad_clip > 0:
                nn.utils.clip_grad_norm_(model.parameters(), max_norm=grad_clip)
            optimizer.step()
        total_loss += float(loss.item())
        steps += 1
    return total_loss / max(steps, 1)

@torch.no_grad()
def evaluate(model: nn.Module, loader: DataLoader, loss_fn: nn.Module, device: torch.device, target_mode: str) -> Dict[str, float]:
    """Evaluate model and return metrics."""
    model.eval()
    total_loss = 0.0
    steps = 0
    preds: List[np.ndarray] = []
    tgts: List[np.ndarray] = []
    for emb, target in loader:
        emb = emb.to(device, non_blocking=True)
        target = target.to(device, non_blocking=True)
        pred = model(emb)
        loss = _compute_loss(pred, target, loss_fn, target_mode)
        total_loss += float(loss.item())
        steps += 1
        preds.append(pred.cpu().numpy())
        tgts.append(target.cpu().numpy())
    pred_np = np.concatenate(preds, axis=0)
    tgt_np = np.concatenate(tgts, axis=0)
    metrics: Dict[str, float] = {'loss': total_loss / max(steps, 1)}
    flat_p = pred_np.reshape(-1)
    flat_t = tgt_np.reshape(-1)
    metrics['mse'] = float(np.mean((flat_t - flat_p) ** 2))
    metrics['mae'] = float(np.mean(np.abs(flat_t - flat_p)))
    pearsons: List[float] = []
    for g in range(pred_np.shape[1]):
        p_g = pred_np[:, g]
        t_g = tgt_np[:, g]
        if np.std(p_g) < 1e-08 or np.std(t_g) < 1e-08:
            continue
        corr = float(np.corrcoef(p_g, t_g)[0, 1])
        if not np.isnan(corr):
            pearsons.append(corr)
    metrics['mean_gene_pearson'] = float(np.mean(pearsons)) if pearsons else 0.0
    return metrics

@torch.no_grad()
def predict_perturbed(model: nn.Module, pert_embeddings: torch.Tensor, target_mode: str, batch_size: int, device: torch.device) -> Dict[str, np.ndarray]:
    """Predict expression for perturbed cells."""
    model.eval()
    n = pert_embeddings.shape[0]
    pred_list: List[np.ndarray] = []
    for start in range(0, n, batch_size):
        end = min(start + batch_size, n)
        batch = pert_embeddings[start:end].to(device, non_blocking=True)
        out = model(batch)
        pred_list.append(out.cpu().numpy())
    result: Dict[str, np.ndarray] = {'pred': np.concatenate(pred_list, axis=0)}
    return result

def get_target_matrix(adata: ad.AnnData, target_mode: str) -> Tuple[Any, str, int]:
    """Return target matrix source, source description, and target gene count."""
    source, source_name = _select_target_source(adata, target_mode)
    n_genes = int(source.shape[1])
    return (source, source_name, n_genes)

def save_outputs(output_dir: Path, args: argparse.Namespace, pred_dict: Dict[str, np.ndarray], train_history: List[Dict[str, Any]], val_metrics: Dict[str, float], best_ckpt: Dict[str, Any], adata_pert: ad.AnnData, pert_indices: np.ndarray, alignment_note: str, emb_shape: Tuple[int, int, int], target_source_name: str, n_bins: Optional[int]) -> None:
    """Save predictions, metrics, config, and metadata artifacts."""
    output_dir.mkdir(parents=True, exist_ok=True)
    pred = pred_dict['pred']
    np.save(output_dir / 'pred_responder_pert.npy', pred)
    torch.save(torch.from_numpy(pred), output_dir / 'pred_responder_pert.pt')
    history_path = output_dir / 'train_history.csv'
    with history_path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(train_history[0].keys()) if train_history else ['epoch'])
        writer.writeheader()
        if train_history:
            writer.writerows(train_history)
    with (output_dir / 'val_metrics.json').open('w') as f:
        json.dump(val_metrics, f, indent=2)
    cfg = vars(args).copy()
    cfg['alignment_note'] = alignment_note
    cfg['best_epoch'] = int(best_ckpt.get('epoch', -1))
    cfg['target_source'] = target_source_name
    with (output_dir / 'train_config.json').open('w') as f:
        json.dump(cfg, f, indent=2)
    with (output_dir / 'best_model.pt').open('wb') as f:
        torch.save(best_ckpt, f)
    with (output_dir / 'used_var_names.txt').open('w') as f:
        for v in adata_pert.var_names.tolist():
            f.write(f'{v}\n')
    obs_df = adata_pert.obs.copy()
    obs_df.insert(0, 'cell_index', pert_indices)
    obs_df.insert(1, 'obs_name', adata_pert.obs_names.tolist())
    obs_df.to_csv(output_dir / 'perturbed_obs_info.csv', index=False)
    shape_summary = {'adata_pert_shape': [int(adata_pert.n_obs), int(adata_pert.n_vars)], 'embedding_shape': list(map(int, emb_shape)), 'prediction_shape': list(map(int, pred.shape)), 'target_mode': args.target_mode, 'target_source': target_source_name, 'n_bins': n_bins, 'alignment_note': alignment_note}
    with (output_dir / 'shape_summary.json').open('w') as f:
        json.dump(shape_summary, f, indent=2)
    with (output_dir / 'shape_summary.txt').open('w') as f:
        for k, v in shape_summary.items():
            f.write(f'{k}: {v}\n')
    try:
        x = pred.astype(np.float32, copy=False)
        adata_pred = ad.AnnData(X=x, obs=adata_pert.obs.copy(), var=adata_pert.var.copy())
        adata_pred.uns['prediction_info'] = {'target_mode': args.target_mode, 'target_source': target_source_name, 'alignment_note': alignment_note, 'n_bins': n_bins}
        adata_pred.write_h5ad(output_dir / 'pred_responder_pert.h5ad')
    except Exception as exc:
        LOGGER.warning('[save] failed to write pred_responder_pert.h5ad: %s', exc)
        np.savez_compressed(output_dir / 'pred_responder_pert_obsvar.npz', X=pred, obs_names=np.asarray(adata_pert.obs_names.tolist(), dtype=object), var_names=np.asarray(adata_pert.var_names.tolist(), dtype=object))

def main(args: argparse.Namespace) -> None:
    """Main training and prediction workflow."""
    configure_logging()
    if args.target_mode != 'data' or args.loss_type != 'huber':
        raise ValueError('Only continuous data targets with Huber loss are supported')
    if not 0.0 < args.train_val_ratio < 1.0:
        raise ValueError('train_val_ratio must be in (0, 1)')
    set_seed(args.seed)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if args.device == 'auto':
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    else:
        device = torch.device(args.device)
    LOGGER.info('[config] device=%s, target_mode=%s', device, args.target_mode)
    adata_path = Path(args.adata_path)
    if not adata_path.exists():
        raise FileNotFoundError(f'AnnData file not found: {args.adata_path}')
    adata = ad.read_h5ad(str(adata_path))
    LOGGER.info('[data] adata shape=%s', adata.shape)
    if 'group' not in adata.obs.columns:
        raise KeyError("adata.obs['group'] is required")
    LOGGER.info('[data] group counts=%s', adata.obs['group'].value_counts().to_dict())
    if Path(args.embedding_path).suffix == '.npy':
        axes_path = getattr(args, 'axes_path', None)
        if not axes_path:
            raise ValueError('NPY embeddings require axes_path from stage 03')
        axes = json.loads(Path(axes_path).read_text())
        if axes.get('leading_cls') is not True or axes.get('obs_names') != adata.obs_names.tolist() or axes.get('var_names') != adata.var_names.tolist():
            raise ValueError('Embedding axes must match AnnData exactly, including leading CLS')
    emb_tensor, emb_gene_names, emb_meta = load_embedding_tensor(args.embedding_path)
    aligned_emb, _, alignment_note = validate_alignment(adata, emb_tensor, emb_gene_names)
    target_source, target_source_name, target_gene_count = get_target_matrix(adata, args.target_mode)
    if target_gene_count != adata.n_vars:
        raise ValueError(f'Target matrix gene dimension does not match adata.n_vars: target_genes={target_gene_count}, adata.n_vars={adata.n_vars}')
    if aligned_emb.shape[1] - 1 != target_gene_count:
        raise ValueError(f'Aligned embedding gene count and target gene count mismatch: emb_genes={aligned_emb.shape[1] - 1}, target_genes={target_gene_count}')
    if sp.issparse(target_source):
        dense_bytes = estimate_dense_bytes(target_source.shape[0], target_source.shape[1], np.float32)
        LOGGER.info('[target] source=%s (sparse), shape=%s, estimated dense float32 size=%.2f GB', target_source_name, target_source.shape, dense_bytes / 1024 ** 3)
    else:
        LOGGER.info('[target] source=%s (dense), shape=%s', target_source_name, target_source.shape)
    pert_mask = adata.obs['group'].values == 'responder_pert'
    nonpert_mask = ~pert_mask
    pert_indices = np.where(pert_mask)[0]
    nonpert_indices = np.where(nonpert_mask)[0]
    if len(pert_indices) == 0:
        raise ValueError("No perturbed cells found where group == 'responder_pert'.")
    if len(nonpert_indices) < 2:
        raise ValueError('Insufficient non-perturbed cells for train/val split.')
    stratify_labels = adata.obs['group'].values[nonpert_indices]
    uniq, cnt = np.unique(stratify_labels, return_counts=True)
    strat_ok = bool(np.all(cnt >= 2))
    if strat_ok:
        try:
            train_idx, val_idx = split_indices(indices=nonpert_indices, labels=stratify_labels, train_ratio=args.train_val_ratio, seed=args.seed, stratified=True)
            LOGGER.info('[split] stratified split succeeded')
        except ValueError as exc:
            LOGGER.warning('[split] stratified split failed (%s); fallback random split', exc)
            train_idx, val_idx = split_indices(indices=nonpert_indices, labels=stratify_labels, train_ratio=args.train_val_ratio, seed=args.seed, stratified=False)
    else:
        LOGGER.warning('[split] stratify not possible due to low group counts=%s; random split', dict(zip(uniq, cnt)))
        train_idx, val_idx = split_indices(indices=nonpert_indices, labels=stratify_labels, train_ratio=args.train_val_ratio, seed=args.seed, stratified=False)
    LOGGER.info('[split] train=%d, val=%d, perturbed=%d', len(train_idx), len(val_idx), len(pert_indices))
    n_bins: Optional[int] = None
    train_emb = aligned_emb[train_idx]
    val_emb = aligned_emb[val_idx]
    pert_emb = aligned_emb[pert_indices]
    train_ds = ExpressionDataset(train_emb, target_source, train_idx, target_gene_count, args.target_mode)
    val_ds = ExpressionDataset(val_emb, target_source, val_idx, target_gene_count, args.target_mode)
    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True, num_workers=args.num_workers, pin_memory=device.type == 'cuda', drop_last=False)
    val_loader = DataLoader(val_ds, batch_size=args.batch_size, shuffle=False, num_workers=args.num_workers, pin_memory=device.type == 'cuda', drop_last=False)
    model = ExpressionPredictor(emb_dim=int(aligned_emb.shape[2]), hidden_dim=args.hidden_dim, target_mode=args.target_mode, n_bins=n_bins, num_layers=args.num_layers, dropout=args.dropout).to(device)
    n_params = sum((p.numel() for p in model.parameters() if p.requires_grad))
    LOGGER.info('[model] trainable parameters=%s', format(n_params, ','))
    loss_fn = nn.HuberLoss()
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    scaler = GradScaler() if args.fp16 and device.type == 'cuda' else None
    best_val_loss = float('inf')
    best_ckpt: Dict[str, Any] = {}
    patience_counter = 0
    train_history: List[Dict[str, Any]] = []
    for epoch in range(1, args.epochs + 1):
        t0 = time.time()
        train_loss = train_one_epoch(model=model, loader=train_loader, optimizer=optimizer, loss_fn=loss_fn, device=device, target_mode=args.target_mode, scaler=scaler, grad_clip=args.grad_clip)
        val_metrics = evaluate(model=model, loader=val_loader, loss_fn=loss_fn, device=device, target_mode=args.target_mode)
        val_loss = float(val_metrics['loss'])
        elapsed = time.time() - t0
        row: Dict[str, Any] = {'epoch': epoch, 'train_loss': train_loss, 'elapsed_sec': elapsed, **val_metrics}
        train_history.append(row)
        LOGGER.info('[epoch %03d/%03d] train_loss=%.6f val_loss=%.6f metrics=%s time=%.1fs', epoch, args.epochs, train_loss, val_loss, {k: v for k, v in val_metrics.items() if k != 'loss'}, elapsed)
        improved = val_loss < best_val_loss - args.early_stop_min_delta
        if improved:
            best_val_loss = val_loss
            patience_counter = 0
            # Preserve live state_dict references for numerical compatibility.
            best_ckpt = {'epoch': epoch, 'model_state_dict': model.state_dict(), 'optimizer_state_dict': optimizer.state_dict(), 'val_loss': val_loss, 'args': vars(args), 'n_bins': n_bins, 'emb_dim': int(aligned_emb.shape[2]), 'alignment_note': alignment_note, 'embedding_meta': emb_meta}
            torch.save(best_ckpt, output_dir / 'best_model.pt')
            LOGGER.info('[checkpoint] updated best model at epoch=%d', epoch)
        else:
            patience_counter += 1
            if patience_counter >= args.patience:
                LOGGER.info('[early_stop] patience reached at epoch=%d', epoch)
                break
    if not best_ckpt:
        raise RuntimeError('Training finished without producing a best checkpoint.')
    model.load_state_dict(best_ckpt['model_state_dict'])
    final_val_metrics = evaluate(model=model, loader=val_loader, loss_fn=loss_fn, device=device, target_mode=args.target_mode)
    pred_dict = predict_perturbed(model=model, pert_embeddings=pert_emb, target_mode=args.target_mode, batch_size=args.batch_size, device=device)
    LOGGER.info('[predict] output prediction shape=%s', pred_dict['pred'].shape)
    adata_pert = adata[pert_indices].copy()
    save_outputs(output_dir=output_dir, args=args, pred_dict=pred_dict, train_history=train_history, val_metrics=final_val_metrics, best_ckpt=best_ckpt, adata_pert=adata_pert, pert_indices=pert_indices, alignment_note=alignment_note, emb_shape=tuple(map(int, aligned_emb.shape)), target_source_name=target_source_name, n_bins=n_bins)
    LOGGER.info('[done] adata shape=%s', adata.shape)
    LOGGER.info('[done] embedding shape=%s', tuple(aligned_emb.shape))
    LOGGER.info('[done] train/val/pert=%d/%d/%d', len(train_idx), len(val_idx), len(pert_indices))
    LOGGER.info('[done] target_mode=%s', args.target_mode)
    LOGGER.info('[done] device=%s', device)
    LOGGER.info('[done] alignment=%s', alignment_note)
    LOGGER.info('[done] output_dir=%s', str(output_dir.resolve()))
