"""Compare paired R/R-prime neighborhoods against R and target-zero NR references.

Input: perturbation/{GENE}_scGPT.h5ad for all requested/manifest genes.
Output: kNN/knn_label_enrichment_{all,filtered}_metrics_*.xlsx/csv.
Preserves the original no-gate, NR gene-zero pairing and empirical statistics.
"""
from argparse import Namespace
from pipeline.common import cli, load_config, output_dir, genes_to_run


def main():
    parser=cli(__doc__)
    parser.add_argument('--genes',nargs='+')
    args=parser.parse_args()
    cfg=load_config(args.config)
    from pipeline import cell_metrics
    if cfg['knn']['nr_condition'] != 'gene_zero':
        raise ValueError('This release preserves only the result-linked NR gene-zero branch')
    out=output_dir(cfg,'kNN')
    selected=genes_to_run(cfg,args.genes)
    # Make discovery independent of leftover result directories, with the original sorted order.
    cell_metrics.discover_gene_dirs=lambda base: [base/g for g in sorted(selected)]
    options=Namespace(base_dir=out,input_dir=output_dir(cfg,'perturbation'),
        min_responder_pert_nonzero=cfg['knn']['min_responder_pert_nonzero'],
        min_nonresponder_gene_zero=cfg['knn']['min_nonresponder_gene_zero_exclusive'],
        target_gene=cfg['knn']['target_gene'], k=cfg['knn']['k'], output_tag=f"k{cfg['knn']['k']}")
    cell_metrics.run(options)


if __name__=='__main__':
    main()
