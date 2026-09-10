"""Set each target gene's responder expression to zero and extract scGPT embeddings.

Input: prepared/IFN_CM_features.h5ad, pretrained configuration and fine-tuned weights.
Output: perturbation/{GENE}_all_concat.h5ad, {GENE}_scGPT.h5ad,
        {GENE}_gene_embeddings.npy (including CLS), {GENE}_axes.json.
--prepare-only creates counterfactuals without GPU inference.
"""
from pipeline.common import cli, load_config, genes_to_run, output_dir, resolve, input_path


def main():
    parser = cli(__doc__)
    parser.add_argument('--genes', nargs='+')
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    cfg = load_config(args.config)
    import scanpy as sc
    from pipeline.perturbation import make_counterfactual, extract_embeddings
    adata = sc.read_h5ad(output_dir(cfg, 'prepared') / 'IFN_CM_features.h5ad')
    out = output_dir(cfg, 'perturbation')
    for gene in genes_to_run(cfg, args.genes):
        combined = make_counterfactual(adata, gene, cfg['metadata']['response'], cfg['metadata']['responder_label'], cfg['perturbation']['target_value'])
        combined.write_h5ad(out / f'{gene}_all_concat.h5ad')
        if not args.prepare_only:
            extract_embeddings(combined, resolve(cfg, cfg['paths']['pretrained_model_dir']), input_path(cfg, 'finetuned_checkpoint'), out / gene, cfg['perturbation'], cfg['metadata'])


if __name__ == '__main__':
    main()
