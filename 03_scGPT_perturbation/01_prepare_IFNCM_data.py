"""Prepare IFN CM inputs from normalized, annotated AnnData.

Input: paths.input_h5ad, pretrained vocab and optionally the archived feature list.
Output: outputs/prepared/IFN_CM_all_genes.h5ad and IFN_CM_features.h5ad.
Purpose: separate the training input from the perturbation feature universe.
The historical training dataset is a separate input; see README.
"""
from pipeline.common import cli, load_config, input_path, output_dir, resolve, require_columns, write_json


def select_features(adata, cfg):
    from scgpt.tokenizer.gene_tokenizer import GeneVocab
    p = cfg['prepare']
    vocab = GeneVocab.from_file(resolve(cfg, cfg['paths']['pretrained_model_dir']) / 'vocab.json')
    if p['feature_mode'] == 'archived':
        genes = input_path(cfg, 'feature_genes').read_text().splitlines()
        if len(genes) != len(set(genes)):
            raise ValueError('Archived feature list contains duplicates')
        missing = [g for g in genes if g not in adata.var_names or g not in vocab]
        if missing:
            raise ValueError(f'Archived features absent from input/vocabulary: {missing[:10]}')
        return adata[:, genes].copy()
    raise ValueError('Only the fixed manuscript feature list is supported')


def main():
    args = cli(__doc__).parse_args()
    cfg = load_config(args.config)
    import scanpy as sc
    m = cfg['metadata']
    adata = sc.read_h5ad(input_path(cfg, 'input_h5ad'))
    require_columns(adata, [m['cell_type'], m['response'], m['batch']])
    keep = adata.obs[m['cell_type']].eq(m['ifn_cm_label'])
    if cfg['prepare']['filter_timepoint']:
        require_columns(adata, [m['timepoint']])
        keep &= adata.obs[m['timepoint']].eq(m['pre_label'])
    adata = adata[keep].copy()
    if not adata.n_obs or adata.obs[m['response']].isna().any():
        raise ValueError('Empty IFN CM subset or missing response labels')
    # Preserve the Seurat export aliases, with current response used only where requested.
    adata.obs['celltype.id'] = adata.obs[m['cell_type']].astype(str).str.replace('+', ' ', regex=False)
    adata.obs['batch.id'] = adata.obs[m['batch']]
    adata.obs['response'] = adata.obs['Response_New'] if 'Response_New' in adata.obs else adata.obs[m['response']]
    adata.var['gene_name'] = adata.var_names
    out = output_dir(cfg, 'prepared')
    adata.write_h5ad(out / 'IFN_CM_all_genes.h5ad')
    selected = select_features(adata, cfg)
    selected.write_h5ad(out / 'IFN_CM_features.h5ad')
    write_json(out / 'preparation_summary.json', {'n_cells':adata.n_obs, 'n_features':selected.n_vars, 'feature_mode':cfg['prepare']['feature_mode'], 'timepoint_filter':cfg['prepare']['filter_timepoint']})


if __name__ == '__main__':
    main()
