"""Construct zero-expression counterfactuals and extract contextual scGPT embeddings.

Inputs: feature-selected AnnData and whole-human configuration/fine-tuned weights.
Outputs: concatenated NR/R/R-prime AnnData, contextual .npy and axis manifest.
The public forward hook exposes the same transformer output as the historical
gene_emb=True modification, without copying or patching the scGPT package.
"""
import json
from pathlib import Path
from pipeline.common import dense, set_seed, write_json


def make_counterfactual(adata, gene, response_column='Response_Final', responder_label='Good', target_value=0):
    import scanpy as sc
    if gene not in adata.var_names:
        raise ValueError(f'Target gene is absent from selected features: {gene}')
    responder = adata[adata.obs[response_column].eq(responder_label)].copy()
    nonresponder = adata[~adata.obs[response_column].eq(responder_label)].copy()
    if not responder.n_obs or not nonresponder.n_obs:
        raise ValueError('Both responder and nonresponder cells are required')
    perturbed = responder.copy()
    idx = perturbed.var_names.get_loc(gene)
    perturbed.obs[f'{gene}_nonzero'] = dense(perturbed.X[:, idx]).ravel() > 0
    # Preserve historical semantics: change X only; do not alter data/raw layers.
    perturbed.X[:, idx] = target_value
    for obj, group in [(nonresponder, 'nonresponder'), (responder, 'responder'), (perturbed, 'responder_pert')]:
        obj.obs['group'] = group
    combined = sc.concat([nonresponder, responder, perturbed], axis=0, join='outer', index_unique='_')
    combined.var['gene_name'] = combined.var_names
    return combined


def extract_embeddings(adata, model_dir, checkpoint, out_prefix, cfg, metadata):
    import numpy as np
    import torch
    from scgpt.model import TransformerModel
    from scgpt.tokenizer.gene_tokenizer import GeneVocab
    from scgpt.preprocess import Preprocessor
    from scgpt.tokenizer import tokenize_and_pad_batch
    set_seed(cfg['seed'])
    if cfg['subset_hvg_on_extraction']:
        raise ValueError('Extraction HVG reselection is disabled in the recorded multi-gene workflow')
    if not cfg['include_zero_gene'] or not cfg['append_cls']:
        raise ValueError('Recorded extraction requires zero genes and a leading CLS token')
    if cfg['expression_source'] != 'X':
        raise ValueError('Recorded extraction reads X')
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    if cfg['fast_transformer'] and device.type != 'cuda':
        raise RuntimeError('Recorded flash-attention inference requires CUDA')
    model_dir = Path(model_dir)
    m = json.loads((model_dir / 'args.json').read_text())
    vocab = GeneVocab.from_file(model_dir / 'vocab.json')
    for special in ['<pad>', '<cls>', '<eoc>']:
        if special not in vocab:
            vocab.append_token(special)
    if any(g not in vocab for g in adata.var_names):
        raise ValueError('Features must be restricted to pretrained vocabulary')
    adata.obs['celltype_id'] = adata.obs[metadata['cell_type']].astype('category').cat.codes.values
    adata.obs['batch_id'] = adata.obs[metadata['batch']].astype('category').cat.codes.values
    adata.var['id_in_vocab'] = 1
    Preprocessor(use_key='X', filter_gene_by_counts=False, filter_cell_by_counts=False,
                 normalize_total=False, subset_hvg=False, binning=cfg['n_bins'],
                 result_binned_key='X_binned')(adata)
    model = TransformerModel(ntoken=len(vocab), d_model=m['embsize'], nhead=m['nheads'],
        d_hid=m['d_hid'], nlayers=m['nlayers'], dropout=m['dropout'],
        pad_token=m['pad_token'], pad_value=m['pad_value'], vocab=vocab,
        input_emb_style=m.get('input_emb_style','continuous'), n_input_bins=m.get('n_bins',51),
        use_fast_transformer=cfg['fast_transformer'], fast_transformer_backend='flash')
    state = torch.load(checkpoint, map_location=device)
    original = model.state_dict()
    matched = {k:v for k,v in state.items() if k in original and v.shape == original[k].shape}
    if not matched:
        raise ValueError('No checkpoint parameters match; verify the checkpoint format')
    # Shape-matched partial loading was also used by the original extraction code.
    original.update(matched)
    model.load_state_dict(original)
    model.to(device).eval()
    ids = np.array(vocab(adata.var_names.tolist()), dtype=int)
    tokenized = tokenize_and_pad_batch(dense(adata.layers['X_binned']), ids, vocab=vocab,
        pad_token=m['pad_token'], pad_value=m['pad_value'], append_cls=True,
        max_len=adata.n_vars+1, include_zero_gene=True)
    # Fail rather than silently labeling shuffled/truncated token positions as genes.
    expected = np.concatenate([[vocab['<cls>']], ids])
    if not np.all(np.asarray(tokenized['genes']) == expected[None, :]):
        raise ValueError('Tokenizer changed feature order; contextual axes would be invalid')
    genes = torch.as_tensor(tokenized['genes'], dtype=torch.long)
    values = torch.as_tensor(tokenized['values'], dtype=torch.float16)
    saved = {}
    def capture(module, inputs, output):
        saved['gene_emb'] = output
    handle = model.transformer_encoder.register_forward_hook(capture)
    contextual = None
    cells = []
    batch_size = cfg['batch_size']
    try:
        with torch.no_grad():
            for start in range(0, adata.n_obs, batch_size):
                batch_g = genes[start:start+batch_size].to(device)
                batch_v = values[start:start+batch_size].to(device)
                with torch.cuda.amp.autocast(enabled=cfg['amp']):
                    result = model(batch_g, batch_v, src_key_padding_mask=batch_g.eq(vocab[m['pad_token']]),
                                   batch_labels=None, CLS=False, CCE=False, MVC=False, ECS=False)
                batch_e = saved.pop('gene_emb').detach().cpu().numpy()
                if contextual is None:
                    contextual = np.lib.format.open_memmap(str(out_prefix)+'_gene_embeddings.npy', mode='w+',
                        dtype=batch_e.dtype, shape=(adata.n_obs, adata.n_vars+1, batch_e.shape[2]))
                contextual[start:start+len(batch_e)] = batch_e
                cells.append(result['cell_emb'].detach().cpu().numpy())
    finally:
        handle.remove()
    contextual.flush()
    adata.obsm['X_scGPT_cell'] = np.concatenate(cells, axis=0)
    adata.write_h5ad(str(out_prefix)+'_scGPT.h5ad')
    write_json(str(out_prefix)+'_axes.json', {'obs_names':adata.obs_names.tolist(), 'var_names':adata.var_names.tolist(), 'leading_cls':True})
    write_json(str(out_prefix)+'_model_loading.json', {'matched_keys':sorted(matched), 'unloaded_keys':sorted(set(original)-set(matched)), 'ignored_checkpoint_keys':sorted(set(state)-set(matched)), 'batch_size':batch_size,'seed':cfg['seed']})
