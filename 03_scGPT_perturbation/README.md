# scGPT in silico perturbation analysis

This directory contains code for simulating gene-expression perturbations in IFN CM cells and comparing responders (R), perturbed responders (R′), and nonresponders (NR).

The workflow uses scGPT to generate contextual gene and cell embeddings. Perturbation effects are evaluated through Euclidean distance changes, gene-level SPC scores, kNN enrichment, expression decoding, AUCell activity, and Gene Ontology Biological Process enrichment.

The workflow supports:

* preparation of annotated expression data and scGPT fine-tuning;
* target-gene perturbation and embedding-based comparisons;
* expression prediction and gene-set activity analysis;
* statistical summaries and visualization.

## Workflow

```text
Annotated single-cell expression data
        ↓
IFN CM preparation and feature selection
        ↓
scGPT fine-tuning
        ↓
Target-gene expression zeroing in responders
        ↓
R / R′ / NR embedding extraction
        ├─ Gene embeddings → SPC and distance changes → Ranking / top genes / GO BP
        ├─ Cell embeddings → kNN enrichment → Ranking
        └─ Expression decoder → Predicted R′ expression → AUCell / random null
```

Existing checkpoints, embeddings, or score tables can enter at the corresponding stage when their input schemas and feature order match.

---

## Installation

Use separate Python/CUDA and R environments on Linux. The Python baseline is Python 3.10.18, scGPT 0.2.4, PyTorch 2.1.0+cu118, and flash-attn 1.0.4. The R baseline is R 4.3.1, Seurat 5.1.0, and AUCell 1.24.0.

| File | Installer | Contents |
|---|---|---|
| `requirements.txt` | Conda | Python, CUDA runtime and build tools with recorded build pins |
| `requirements-pip.txt` | pip | Numerical libraries and their pinned runtime dependencies |
| `requirements-r.txt` | Conda | R and packages for AUCell and GO BP |

Create the Python environment:

```bash
conda create -n scgpt-perturbation --override-channels \
  -c conda-forge --file requirements.txt
conda activate scgpt-perturbation

python -m pip install -r requirements-pip.txt
python -m pip install --no-deps scgpt==0.2.4

export CUDA_HOME="$CONDA_PREFIX"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
MAX_JOBS=4 python -m pip install --no-build-isolation --no-deps flash-attn==1.0.4
```

`requirements.txt` is a Conda specification, not a pip requirements file. scGPT is installed with `--no-deps` after the runtime dependencies because its package metadata also requests unrelated scVI/JAX workflows. This environment targets the scripts in this repository. The library-path setting selects the Conda C++ runtime required by the compiled attention extension. The recorded CUDA runtime is 11.8 and development toolkit package is 11.7; flash-attention compilation also depends on the host NVIDIA driver and GPU. See the [scGPT installation instructions](https://github.com/bowang-lab/scGPT#installation) for platform requirements.

Create the R environment and connect it to the active Python environment:

```bash
conda create -n scgpt-perturbation-r --override-channels \
  -c conda-forge -c bioconda --file requirements-r.txt

# Run while the scgpt-perturbation Python environment is active.
export RETICULATE_PYTHON="$(command -v python)"
```

The R requirements pin installed package versions; they are not a complete transitive lockfile. These specifications were checked against installed environments, but a clean installation and GPU rebuild have not been tested. The AUCell script records R session information in its output directory.

Obtain the whole-human pretrained bundle from the [scGPT model zoo](https://github.com/bowang-lab/scGPT#pretrained-scgpt-model-zoo). Place matching `args.json`, `vocab.json`, and `best_model.pt` files in the directory configured by `paths.pretrained_model_dir`.

---

## Inputs and configuration

Copy the example configuration and edit the paths and metadata mapping:

```bash
cp config.example.yaml config.yaml
```

Main inputs are:

| Input | Description |
|---|---|
| Expression data | Normalized AnnData H5AD with unique cell identifiers and gene symbols |
| Training data | IFN CM H5AD with `celltype.id` and `batch.id`, when using a separate fine-tuning input |
| Pretrained model | Matching scGPT `args.json`, `vocab.json`, and `best_model.pt` |
| Feature list | Ordered gene symbols used for perturbation embedding extraction |
| Perturbation manifest | TSV with `gene` and `group_label` columns |
| AUCell reference | Annotated Seurat RDS containing normalized RNA expression |
| Decoder predictions | R′ expression H5AD produced by stage 06 |
| Distance-change scores | Workbook configured through `paths.distance_change_scores` |
| GO BP results | Cached EnrichR workbook for stage 09, or results obtained with `--query` |

Map the input metadata columns and category labels under `metadata`. Required information includes cell type, response group, batch, and timepoint. The configured responder label identifies R; other nonmissing response labels identify NR. AUCell selects the configured pre-treatment IFN CM population.

Preparation reads normalized expression from `X` without additional normalization or log transformation. Preserve `layers['data']` for decoder training; the decoder uses `X` when that layer is absent. Set `finetune.use_prepared_input: true` to train on the stage-01 input, or supply a separate training dataset through `paths.finetune_input_h5ad`.

Stage 03 generates group labels, the pre-perturbation `{GENE}_nonzero` flag, cell embeddings in `obsm['X_scGPT_cell']`, and binned expression in `layers['X_binned']`.

The distance-change workbook must contain an `All_Results` sheet with `pre_euclidean_distance_mean` and `post_euclidean_distance_mean`. Configure this input independently of the SPC ranking tables.

All relative paths resolve against the configuration file's directory. Analysis parameters and input paths are defined in [config.example.yaml](config.example.yaml).

---

## Fine-tuning and simulated perturbation

Fine-tuning uses a shuffled 90%/10% cell split, binned expression, masked expression reconstruction, MVC/GEPC reconstruction, explicit-zero objectives, ECS and batch-adversarial loss. Model and optimizer settings are defined in `finetune.hparams`.

For each target, responder cells are copied and the target's expression in `X` is set to zero. The pre-perturbation positive-expression flag is retained. Other expression layers are preserved. NR, R and R′ are concatenated, binned and passed through scGPT with zero-expression genes and a leading CLS token. This simulates an input-expression change; it does not estimate a causal intervention effect.

Contextual embeddings have shape `cells × (genes + CLS) × embedding dimension`. NPY exports include an axis manifest to validate cell and gene order. Cell embeddings are stored in the accompanying H5AD.

---

## SPC and kNN enrichment

The score implementation is named `SRC` in the entry point and stored as `SPC` in output columns. Gene-level SPC compares target-positive R and corresponding R′ with all NR. After pooling contextual embeddings across cells, define:

```text
SPC = 2 × (d_pre − d_post) / (d_pre + d_post)
```

Here, `d_pre` is the R-to-NR distance and `d_post` is the R′-to-NR distance for a gene. A zero denominator produces a score of zero. Positive values indicate movement toward NR. The perturbation target and CLS token are excluded from scored genes. The ranking uses the average of the 20 largest mean-pooled Euclidean SPC values. The distance-change panel separately plots `d_pre − d_post` for every scored gene as a 100-bin histogram with a logarithmic count axis; it does not use SPC normalization.

kNN enrichment uses paired target-positive R/R′ and NR cells with zero target expression. Neighbors are fitted on R and eligible NR. For each pair, the score is `NR_fraction(R′) − NR_fraction(R)`, excluding R's own index from its reference query. The primary ranking uses the mean Euclidean shift with `k=100`.

Perturbation rankings compare a target with the other eligible perturbations using `(1 + count(background ≥ target)) / (1 + background size)`. Cell-count thresholds are configured separately for SPC and kNN.

---

## Expression decoding and AUCell

The decoder trains on unperturbed R and NR, using a group-stratified 90%/10% split when possible. Each gene's input concatenates its contextual embedding, the cell CLS embedding, and their elementwise product. A shared MLP uses LayerNorm, two hidden layers with GELU and dropout, and a scalar output. The configured continuous-expression model uses Huber loss and AdamW. Predictions are not clipped or converted to counts.

Stage 06 writes `pred_responder_pert.h5ad`, with R′ cells × genes in `X` and corresponding observation/variable metadata. It also saves NPY/PT predictions, a checkpoint, training history and validation metrics. These predictions are needed for AUCell; SPC and kNN operate directly on embeddings.

AUCell ranks original and predicted expression separately over their common gene universe. For paired cells, `delta = AUC_predicted − AUC_original`. The analysis reports negative, zero and positive delta proportions. Equally sized random gene sets provide a null distribution of trimmed-mean deltas; its lower-tail p-value is `mean(null ≤ observed)`.

The analysis uses **no detection-frequency filtering**, **1,000 random gene sets**, a 10% trimmed mean and null sampling seed 1. The random pool is the common gene universe excluding SIGLEC1. Both original and predicted matrices are ranked over their common genes.

---

## Perturbation background and GO BP

The perturbation manifest defines the gene universe used for background ranking. Group labels are retained as annotations; ranking calculations use all eligible perturbations without grouping by housekeeping status. Gene-list sampling is performed outside this workflow.

GO BP enrichment uses the SPC top-20 genes for each perturbation and the `GO_Biological_Process_2021` database. Terms are ordered by adjusted p-value and overlap fraction, and the top ten are plotted. The adjusted-p-value reference line is not a term-exclusion threshold.

Stage 09 reads cached EnrichR results by default. Use `--query` to obtain results from the live service.

---

## Main parameters

Important parameters are defined in `config.example.yaml`.

| Parameter | Default | Description |
|---|---:|---|
| `finetune.hparams.n_hvg` | `3000` | Requested training HVGs |
| `finetune.hparams.n_bins` | `51` | Expression bins |
| `finetune.hparams.epochs` | `50` | Fine-tuning epochs |
| `finetune.hparams.lr` | `0.0001` | Fine-tuning learning rate |
| `finetune.hparams.seed` | `0` | Fine-tuning random seed |
| `perturbation.target_value` | `0` | Simulated target expression |
| `perturbation.seed` | `42` | Embedding extraction seed |
| `src.primary_top_n` | `20` | Genes used for the SPC ranking summary |
| `src.min_original_nonzero` | `207` | Minimum original-positive R′ cells for SPC ranking |
| `knn.k` | `100` | Nearest neighbors |
| `knn.min_responder_pert_nonzero` | `135` | Minimum original-positive R′ cells for kNN ranking |
| `knn.min_nonresponder_gene_zero_exclusive` | `117` | Strict lower bound on target-zero NR cells |
| `decoder.num_layers` | `2` | Hidden decoder layers |
| `decoder.hidden_dim` | `512` | Units per hidden layer |
| `decoder.epochs` | `200` | Maximum decoder training epochs |
| `decoder.seed` | `42` | Decoder random seed |
| `aucell.filter` | `nofilter` | Gene detection-frequency filtering |
| `aucell.n_iterations` | `1000` | Random gene-set iterations |
| `aucell.trim` | `0.1` | Trimming fraction per tail |
| `aucell.seed_base` | `1` | Random gene-set sampling seed |
| `gobp.top_terms` | `10` | GO BP terms displayed |

AUCell ranking tie-breaking is unseeded when `aucell.ranking_seed` is null. Sampling seed 1 alone therefore does not fix all random-number behavior. Model results can also depend on package versions, GPU topology, batch size, and input order.

---

## Model loading and checkpoint behavior

Embedding extraction loads shape-compatible checkpoint parameters and records unmatched keys. Use matching model configurations and vocabularies when supplying existing weights.

The decoder retains live checkpoint tensor references during training. Its final saved weights can therefore correspond to a later epoch than the stored best-epoch label. This behavior is preserved in the implementation.

---

## Scripts

| Script | Purpose |
|---|---|
| `01_prepare_IFNCM_data.py` | Prepare IFN CM expression data and the fixed feature universe |
| `02_finetune_scGPT.py` | Fine-tune the pretrained model |
| `03_run_perturbation.py` | Create perturbations and extract embeddings |
| `04_SRC_analysis.py` | Calculate SPC and top-20 gene summaries |
| `05_kNN_enrichment.py` | Calculate Euclidean kNN enrichment |
| `06_expression_decoder.py` | Train the continuous-expression decoder and predict R′ expression |
| `07_AUCell_analysis.R` | Calculate AUCell delta signs and the random-null distribution |
| `08_plot_results.py` | Plot SPC ranks, kNN ranks, top-20 genes, and distance changes |
| `09_GO_BP_enrichment.R` | Calculate and plot GO BP enrichment |
| `pipeline/` | Shared training, analysis, plotting, and configuration functions |
| `resources/` | Feature and perturbation gene lists |
| `tests/test_analysis.py` | Synthetic numerical regression tests |

---

## Running the main workflow

```bash
cp config.example.yaml config.yaml

# Edit config.yaml before running the analysis.

python 01_prepare_IFNCM_data.py --config config.yaml

# Set the process count to the available GPUs.
torchrun --standalone --nproc_per_node=1 02_finetune_scGPT.py --config config.yaml
```

Set `paths.finetuned_checkpoint` to the checkpoint selected for embedding extraction, then continue:

```bash
python 03_run_perturbation.py --config config.yaml
python 04_SRC_analysis.py --config config.yaml
python 05_kNN_enrichment.py --config config.yaml
python 06_expression_decoder.py --config config.yaml --genes SIGLEC1
conda run --no-capture-output -n scgpt-perturbation-r Rscript 07_AUCell_analysis.R config.yaml
python 08_plot_results.py --config config.yaml
conda run --no-capture-output -n scgpt-perturbation-r Rscript 09_GO_BP_enrichment.R config.yaml
```

Stages 03–05 use the full perturbation manifest by default. Passing `--genes` restricts execution and the background available for ranking. The plotting routines highlight SIGLEC1, SERPINB2, and THBD; the AUCell implementation is specific to SIGLEC1.

Available execution options include:

| Option | Purpose |
|---|---|
| Stage 03: `--prepare-only` | Create counterfactual inputs without embedding extraction |
| Stage 04: `--tables-only` | Summarize existing gene-score workbooks |
| Stage 06: input/output path arguments | Use existing embeddings for a single target |
| Stage 09: `--query` | Request GO BP results from EnrichR |

Run Python entry points with `--help` for their arguments. Use a separate output directory for each configuration.

Synthetic numerical checks can be run with:

```bash
python -m unittest discover -s tests -v
```

---

## Outputs

The workflow produces:

* prepared expression data, embeddings, and fine-tuned checkpoints;
* SPC and kNN summary tables and ranked-distribution figures;
* top-20 gene barplots and the gene-distance-change histogram;
* predicted R′ expression and decoder validation metrics;
* AUCell delta-sign proportions and the nofilter/n1000 random-null figure;
* GO BP enrichment tables and figures.

Tables are written as CSV/XLSX, figures as PDF/PNG, and expression matrices as H5AD. Contextual NPY embeddings include JSON axis metadata.

Input datasets, pretrained weights, trained checkpoints, and generated outputs are not included in this repository. Local configurations, data, model files, and outputs are excluded through `.gitignore`.
