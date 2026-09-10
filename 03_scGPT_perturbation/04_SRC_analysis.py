"""Compute gene-wise SRC/SPC and rank perturbations by their top-N gene scores.

Input: perturbation/{GENE}_scGPT.h5ad and contextual embeddings (NPY or original PT).
Output: SRC/{GENE}/nonresponder_similarity_analysis.xlsx, pooled score matrix,
        top-N selected-gene tables and the >=207 original-positive-cell ranking table.
Use --tables-only to summarize previously computed per-gene XLSX files.
"""
import json
from pipeline.common import cli, load_config, genes_to_run, output_dir, input_path


def main():
    parser = cli(__doc__)
    parser.add_argument('--genes', nargs='+')
    parser.add_argument('--tables-only', action='store_true')
    args = parser.parse_args()
    cfg = load_config(args.config)
    if (cfg['src']['primary_top_n'] != 20 or cfg['src']['primary_pooling'] != 'mean'
            or cfg['src']['primary_distance'] != 'euclidean' or not cfg['src']['exclude_target_gene']
            or cfg['src']['top_n'] != [20]):
        raise ValueError('The preserved Figure 5 definition requires top-20 mean Euclidean SPC, excluding target')
    import numpy as np
    import pandas as pd
    import anndata as ad
    from pipeline.src_metrics import analyze_gene, summarize_table, group_indices
    out, inp = output_dir(cfg, 'SRC'), output_dir(cfg, 'perturbation')
    rows, selected = [], {k:[] for k in cfg['src']['top_n']}
    labels = pd.read_csv(input_path(cfg,'gene_manifest'),sep='\t').set_index('gene').group_label.to_dict()
    for gene in genes_to_run(cfg, args.genes):
        directory = out / gene
        directory.mkdir(exist_ok=True)
        if args.tables_only:
            table = pd.read_excel(directory / 'nonresponder_similarity_analysis.xlsx',sheet_name=0)
            adata = ad.read_h5ad(inp / f'{gene}_scGPT.h5ad')
            n = len(group_indices(adata,gene)[1])
        else:
            adata = ad.read_h5ad(inp / f'{gene}_scGPT.h5ad')
            npy = inp / f'{gene}_gene_embeddings.npy'
            axes = inp / f'{gene}_axes.json'
            if axes.exists():
                names=json.loads(axes.read_text())
                if names['obs_names'] != adata.obs_names.tolist() or names['var_names'] != adata.var_names.tolist():
                    raise ValueError('Embedding manifest does not match AnnData order')
            if npy.exists():
                if not axes.exists():
                    raise ValueError('NPY embeddings require an axis manifest')
                emb = np.load(npy,mmap_mode='r')
            else:
                import torch
                # Original tensor-only artifacts have no identifiers; preserve their documented order assumption.
                emb = torch.load(inp / f'{gene}_gene_embeddings.pt',map_location='cpu',mmap=True).numpy()
            table,n = analyze_gene(adata,emb,gene)
            table.to_excel(directory/'nonresponder_similarity_analysis.xlsx',sheet_name='All_Results',index=False)
            del emb
        row, tops = summarize_table(table,gene,cfg['src']['top_n'])
        row['group_label'] = labels.get(gene,gene)
        row['n_original_nonzero'] = n
        rows.append(row)
        for k in selected:
            selected[k].append(tops[k])
    matrix = pd.DataFrame(rows).sort_values('gene')
    matrix.to_csv(out/'spc_mean_matrix.csv',index=False)
    matrix.to_excel(out/'spc_mean_matrix.xlsx',index=False)
    filtered = matrix[matrix.n_original_nonzero >= cfg['src']['min_original_nonzero']]
    filtered.to_excel(out/f"spc_mean_matrix_true_count_ge{cfg['src']['min_original_nonzero']}.xlsx",index=False)
    for k, records in selected.items():
        df=pd.DataFrame(records).sort_values('perturb_gene')
        df.to_csv(out/f'spc_euclidean_top{k}_selected_genes.csv',index=False)
        df.to_excel(out/f'spc_euclidean_top{k}_selected_genes.xlsx',index=False)


if __name__ == '__main__':
    main()
