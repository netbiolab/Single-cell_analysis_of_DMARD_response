"""Shared configuration, paths, random seeds and input validation."""
from pathlib import Path
import argparse
import csv
import json


def load_config(path):
    import yaml
    path = Path(path).resolve()
    cfg = yaml.safe_load(path.read_text())
    cfg['_config_dir'] = path.parent
    return cfg


def cli(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--config', required=True)
    return parser


def resolve(cfg, value):
    path = Path(value).expanduser()
    return path if path.is_absolute() else cfg['_config_dir'] / path


def input_path(cfg, key):
    path = resolve(cfg, cfg['paths'][key])
    if not path.is_file():
        raise FileNotFoundError(f'Missing paths.{key}: {path}')
    return path


def output_dir(cfg, stage):
    path = resolve(cfg, cfg['paths']['output_dir']) / stage
    path.mkdir(parents=True, exist_ok=True)
    return path


def set_seed(seed):
    # Exactly the RNGs set by the historical private helper.
    import random
    import numpy as np
    import torch
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def genes_to_run(cfg, requested=None):
    if requested:
        genes = requested
    elif cfg['perturbation']['all_manifest_genes']:
        with input_path(cfg, 'gene_manifest').open() as f:
            genes = [r['gene'] for r in csv.DictReader(f, delimiter='\t')]
    else:
        genes = cfg['perturbation']['genes']
    if len(set(genes)) != len(genes):
        raise ValueError('Duplicate perturbation genes')
    for gene in genes:
        if not gene or '/' in gene or '\\' in gene or gene in {'.', '..'}:
            raise ValueError('Invalid gene name')
    return genes


def dense(x):
    import numpy as np
    return x.toarray() if hasattr(x, 'toarray') else np.asarray(x)


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, default=str)+'\n')


def require_columns(adata, columns):
    missing = set(columns) - set(adata.obs.columns)
    if missing:
        raise ValueError(f'Missing obs columns: {sorted(missing)}')
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError('Cell and gene names must be unique before concatenation')
