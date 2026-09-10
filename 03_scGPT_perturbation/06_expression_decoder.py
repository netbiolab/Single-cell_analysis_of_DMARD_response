"""Decode R-prime contextual embeddings into expression for AUCell.

INPUT: stage-03 {GENE}_scGPT.h5ad, CLS + gene embeddings (NPY + axes JSON
or trusted historical PT), and YAML decoder settings.
OUTPUT: outputs/decoder/{GENE}/pred_responder_pert.h5ad (cells x genes),
NPY/PT copies, checkpoint, train history, validation metrics and gene/cell axes.
"""
import argparse
from pipeline.common import cli, load_config, resolve, genes_to_run


def main():
    parser = cli(__doc__)
    parser.add_argument('--genes', nargs='+')
    for name in ('adata-path', 'embedding-path', 'axes-path', 'output-dir'):
        parser.add_argument('--'+name, help='Single-gene path override; relative to config')
    opts = parser.parse_args()
    cfg = load_config(opts.config)
    settings = cfg['decoder']
    genes = genes_to_run(cfg, opts.genes or settings['genes'])
    if len(genes) != 1 and any((opts.adata_path, opts.embedding_path, opts.axes_path, opts.output_dir)):
        parser.error('Path overrides require exactly one --genes value')
    keys = ('target_mode batch_size epochs lr weight_decay hidden_dim dropout '
            'num_layers patience seed num_workers device loss_type '
            'fp16 train_val_ratio grad_clip early_stop_min_delta').split()
    root = resolve(cfg, cfg['paths']['output_dir'])
    jobs = []
    for gene in genes:
        stage = root / 'perturbation'
        emb = resolve(cfg, opts.embedding_path) if opts.embedding_path else stage/(gene+'_gene_embeddings.npy')
        if not opts.embedding_path and not emb.is_file():
            emb = emb.with_suffix('.pt')
        adata = resolve(cfg, opts.adata_path) if opts.adata_path else stage/(gene+'_scGPT.h5ad')
        axes = resolve(cfg, opts.axes_path) if opts.axes_path else stage/(gene+'_axes.json')
        for path in [adata, emb] + ([axes] if emb.suffix == '.npy' else []):
            if not path.is_file():
                raise FileNotFoundError(path)
        destination = resolve(cfg, opts.output_dir) if opts.output_dir else root/'decoder'/gene
        jobs.append(argparse.Namespace(**{k: settings[k] for k in keys},
                    adata_path=str(adata), embedding_path=str(emb),
                    axes_path=str(axes) if emb.suffix == '.npy' else None,
                    output_dir=str(destination)))
    from pipeline.expression_decoder import main as train_predict
    for args in jobs:
        train_predict(args)


if __name__ == '__main__':
    main()
