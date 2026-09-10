"""Gene-wise SRC (historical column prefix: SPC) and exact original distance kernels."""
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances
from pipeline.common import dense


def apply_pooling(embeddings):
    return np.mean(embeddings, axis=0)


def calculate_gene_distances_vs_target(emb1, target_emb, gene_indices):
    return np.array([euclidean_distances(emb1[i].reshape(1, -1),
                    target_emb[i].reshape(1, -1))[0, 0] for i in gene_indices])


def calculate_spc(pre, post):
    # The source docstring gave the opposite sign; its executable expression is retained.
    denominator = post + pre
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(denominator != 0, 2 * (pre - post) / denominator, 0)


def group_indices(adata, gene):
    groups = adata.obs['group'].astype(str).to_numpy()
    expression = dense(adata[:, gene].X).ravel()
    flag = f'{gene}_nonzero'
    if flag not in adata.obs:
        raise ValueError(f'Missing pre-perturbation flag {flag}')
    # Matches the tensor template: R gated by X>0, R-prime by the saved flag, all NR.
    r = np.flatnonzero((groups == 'responder') & (expression > 0))
    rp = np.flatnonzero((groups == 'responder_pert') & adata.obs[flag].eq(True).to_numpy())
    nr = np.flatnonzero(groups == 'nonresponder')
    if min(len(r), len(rp), len(nr)) == 0:
        raise ValueError(f'Empty SRC group for {gene}')
    return r, rp, nr


def analyze_gene(adata, embeddings, gene):
    if embeddings.ndim != 3 or embeddings.shape[0] != adata.n_obs:
        raise ValueError('Contextual embedding cell axis does not match AnnData')
    if embeddings.shape[1] == adata.n_vars+1:
        embeddings = embeddings[:, 1:, :]
    elif embeddings.shape[1] != adata.n_vars:
        raise ValueError('Contextual gene axis does not match AnnData')
    r, rp, nr = group_indices(adata, gene)
    # Materialization and dtype of pooling intentionally match the original template.
    by_group = [embeddings[idx] for idx in (r, rp, nr)]
    gene_indices = [i for i, name in enumerate(adata.var_names) if name != gene]
    result = {'gene_name':[adata.var_names[i] for i in gene_indices]}
    a, b, c = [apply_pooling(e) for e in by_group]
    pre = calculate_gene_distances_vs_target(a, c, gene_indices)
    post = calculate_gene_distances_vs_target(b, c, gene_indices)
    result['pre_euclidean_distance_mean'] = pre
    result['post_euclidean_distance_mean'] = post
    result['change_euclidean_mean'] = pre-post
    result['spc_euclidean_mean'] = calculate_spc(pre, post)
    return pd.DataFrame(result), len(rp)


def summarize_table(table, gene, top_ns):
    row = {'gene':gene}
    for metric in ['spc_euclidean_mean']:
        row[metric] = table[metric].mean()
    selected = {}
    for k in top_ns:
        top = table.nlargest(k, 'spc_euclidean_mean')[['gene_name','spc_euclidean_mean']].reset_index(drop=True)
        mean = top.spc_euclidean_mean.mean()
        row[f'top{k}_mean_spc_euclidean_mean'] = mean
        entry = {'perturb_gene':gene, f'top{k}_mean_spc_euclidean_mean':mean, 'selected_count':len(top)}
        for i in range(k):
            entry[f'top{i+1}_gene'] = top.iloc[i].gene_name if i < len(top) else None
            entry[f'top{i+1}_spc_euclidean_mean'] = top.iloc[i].spc_euclidean_mean if i < len(top) else None
        selected[k] = entry
    return row, selected
