# Responder Score analysis

This directory contains code for calculating IFNα−IFNγ reference-based Responder Scores and comparing responders (R) with nonresponders (NR).

Scores are calculated from mean expression profiles for each patient within each fine cell type. The workflow includes reference projection, pooled score-distribution comparison and a cell-type-resolved rank-biserial analysis with exact synchronized patient-label permutation.

## Workflow

```text
Annotated single-cell expression data
        ↓
Patient × fine-cell-type mean expression
        ↓
IFNα and IFNγ reference projection
        ↓
Responder Score calculation
        ↓
R / NR distribution and patient-level comparisons
        ↓
Statistical summaries and visualization
```

An existing score table can enter directly at the statistical-analysis stage after input validation.

---

## Installation

Create the analysis environment using [requirements.txt](requirements.txt):

```bash
conda create -n responder_score --override-channels \
  -c conda-forge --file requirements.txt

conda activate responder_score
```

The requirements file uses conda package specifications for Linux (`linux-64`); it is not a pip requirements file. Package versions and builds are pinned for the analysis dependencies and their resolved dependencies.

| Package | Version |
| --- | --- |
| R | 4.3.1 |
| Seurat | 5.1.0 |
| SeuratObject | 5.0.2 |
| Matrix | 1.6-5 |
| dplyr | 1.1.4 |
| yaml | 2.3.7 |
| ggplot2 | 3.5.1 |

Statistical analysis of an existing score table requires base R and yaml; plotting additionally requires ggplot2 and Cairo PDF support. Each stage writes a package-version record. Package provenance and compatibility checks are documented in [ENVIRONMENT.md](ENVIRONMENT.md).

---

## Inputs and configuration

Copy the example configuration and edit the paths and metadata mapping:

```bash
cp config.example.yaml config.yaml
```

| Input | Description |
| --- | --- |
| Expression data | Annotated Seurat RDS containing normalized expression in `RNA/data` |
| Metadata | Patient identifier, fine cell type, and response group for each cell |
| Reference weights | Included `resources/centroids_top7000.tsv`, with `axis`, `rank`, `gene`, and `weight` columns |
| Cytokine DEG table | Optional CSV with `cytokine`, `adj_p_value`, `gene`, and `log_fc`, for rebuilding the references |

Input locations are specified by `paths.seurat`, `paths.centroids`, and `paths.cytokine_degs`; results are written to `paths.output`. Relative paths resolve from the configuration file's directory.

Map the input metadata through `columns.patient`, `columns.celltype`, and `columns.response`. Set `labels.R` and `labels.NR` to the corresponding response values. Each patient must have a consistent response label.

To analyze existing scores, place `projection_scores.csv` in the configured output directory. It must contain `donor`, `scvi_fine_CT`, `response`, and `Responder_score_cosine`, with one finite score per patient × cell type. Projection scatter plots also require `IFNa_cosine`, `IFNg_cosine`, `IFNa_correlation`, and `IFNg_correlation`.

The exact-test implementation includes cohort-size, group-size, cell-type-count, and permutation-evaluability checks. It is configured for the accompanying analysis design; applying it to another cohort requires reviewing these configuration values and the cell-type availability assertions in `R/patient_robustness.R`.

---

## Expression aggregation and reference construction

For patient `p`, fine cell type `t`, and gene `g`, calculate:

```text
x[p,t,g] = mean(RNA/data[g, cells from patient p and cell type t])
```

Groups containing fewer than the configured minimum number of cells are omitted. The supplied normalized expression is averaged directly, without additional normalization or gene-wise standardization. The score is calculated from the resulting patient × cell-type profile.

The default `centroid.source: frozen` loads the included ordered reference weights. To construct references from a cytokine DEG table, set `centroid.source: deg_table`.

For each reference cytokine (`IFN-alpha1` and `IFN-gamma`):

1. Retain records with adjusted p-value below 0.05.
2. Average `log_fc` by gene across retained records.
3. Sort genes by decreasing mean log fold change.
4. Select the first 7,000 genes and use their mean log fold changes as weights.

Weights retain their signs; no additional positive-log-fold-change filter is applied. Reference genes are selected from the cytokine DEG table, without selection on patient response labels.

---

## Responder Score

For each reference axis `k`, use the intersection `G_k` of expression genes and reference genes:

```text
cos_k(x) = sum(x[g] * w_k[g]) /
           (sqrt(sum(x[g]^2)) * sqrt(sum(w_k[g]^2)))
```

All sums are over `G_k`. Gene intersections and vector norms are calculated independently for the two axes.

```text
Responder_score_cosine = cos_IFNa(x) - cos_IFNg(x)
```

The auxiliary `Responder_score_correlation` is the difference between the two Pearson expression/reference correlations. Positive scores are classified as IFNα-like; zero or negative scores are classified as IFNγ-like.

Fewer than the configured minimum number of shared genes returns a missing score. Statistical analysis requires finite scores and stops if this requirement is not met.

---

## Statistical analysis

### Pooled score-distribution comparison

A two-sided Kolmogorov–Smirnov test compares R and NR cosine-score distributions across patient × cell-type profiles. The implementation uses `ks.test(..., exact = NULL)`, allowing R to select its default exact or asymptotic calculation.

Profiles from the same patient share the patient as their sampling unit. The pooled test describes score distributions; the following analysis uses patients as the label-permutation unit.

### Cell-type-resolved rank-biserial comparison

Within each fine cell type, rank the available patient scores using average ranks for ties. Calculate:

```text
U_R  = sum(ranks of R patients) - n_R * (n_R + 1) / 2
r_rb = 2 * U_R / (n_R * n_NR) - 1
```

The rank-biserial correlation ranges from −1 to 1. Positive values indicate higher scores in responders. The corresponding pairwise probability, with half credit for ties, is `(r_rb + 1) / 2`.

Aggregate cell-type effects using a 10% trimmed mean. Sort the effects, remove `floor(0.10 * number of evaluable types)` from each tail, and equally average the remaining effects.

The exact synchronized permutation test proceeds as follows:

1. Enumerate all patient-label assignments preserving the R/NR group sizes.
2. Apply each patient's assigned label to every cell type, keeping scores and missing-type structure fixed.
3. Recalculate the within-type rank-biserial effects. Omit a type for that assignment if either response group is absent.
4. Recalculate the trimmed mean, including a new sorting and trimming step for each assignment.
5. Calculate the one-sided upper-tail fraction using:

```text
tolerance = sqrt(.Machine$double.eps) * max(1, abs(T_observed))
T_permuted >= T_observed - tolerance
```

All assignments, including the observed assignment, remain in the denominator. The excluded cell types may differ between assignments. This comparison is controlled by `statistics.rank_biserial_permutation`.

The exact test uses complete enumeration, without a +1 correction, bootstrap, or multiple-testing adjustment. No random sampling is performed, so no random seed is required.

---

## Main parameters

Analysis and plotting parameters are defined in [config.example.yaml](config.example.yaml).

| Parameter | Default | Description |
| --- | --- | --- |
| `centroid.source` | `frozen` | Load included reference weights |
| `centroid.adjusted_p_cutoff` | `0.05` | Reference DEG significance threshold |
| `centroid.top_n` | `7000` | Genes retained per reference axis |
| `centroid.minimum_cells` | `10` | Minimum cells per patient × cell-type profile |
| `centroid.minimum_shared_genes` | `10` | Minimum shared genes per projection axis |
| `statistics.permutations` | `6435` | Expected number of completely enumerated assignments |
| `rank_biserial.trim_fraction` | `0.1` | Fraction removed from each tail |
| `rank_biserial.alternative` | `greater` | Responder-positive upper-tail test |
| `statistics.seed` | `null` | No RNG used |
| `plot.bins` | `30` | Score-distribution histogram bins |

Cohort and evaluability expectations are recorded under `statistics.expected_*`, `rank_biserial.expected_*`, and `rank_biserial.permutations_with_*`. Estimator descriptors such as the aggregation method and test direction document the implemented methods; they are not general-purpose switches for alternative analyses.

---

## Scripts

| Script | Purpose |
| --- | --- |
| `01_prepare_data.R` | Validate metadata and prepare patient × cell-type mean expression |
| `02_calculate_responder_score.R` | Load or build references and calculate cosine/Pearson scores |
| `03_statistical_analysis.R` | Run pooled KS and trimmed rank-biserial comparisons |
| `04_plot_results.R` | Generate score-distribution and reference-projection figures |
| `R/common.R` | Shared configuration, reference, projection, and aggregation functions |
| `R/patient_robustness.R` | Patient-matrix validation and the synchronized rank-biserial exact test |
| `tests/test_numerics.R` | Synthetic numerical checks |

---

## Running the workflow

Run the following commands from this directory:

```bash
cp config.example.yaml config.yaml

# Edit config.yaml before running the analysis.

Rscript 01_prepare_data.R --config config.yaml
Rscript 02_calculate_responder_score.R --config config.yaml
Rscript 03_statistical_analysis.R --config config.yaml
Rscript 04_plot_results.R --config config.yaml
```

For an existing score table, start at stage 03. Stage 04 requires the pooled KS output, so keep `statistics.group_comparison: true` when generating the figures. Use a separate output directory for each configuration.

Synthetic numerical checks can be run with:

```bash
Rscript tests/test_numerics.R
```

---

## Outputs

| Output | Description |
| --- | --- |
| `prepared.rds` | Mean-expression matrix and patient × cell-type metadata |
| `projection_scores.csv` | Reference projections, scores, and classifications |
| `selected_genes.csv` | Reference genes, weights, ordering, and presence in expression data |
| `statistical_tests.csv`, `group_summary.csv` | Pooled KS result and response-group summaries |
| `rank_biserial_celltype_effects.csv` | Within-type effects, pairwise probabilities, and observed trim status |
| `rank_biserial_trimmed_mean_exact_summary.csv` | Trimmed effect and one-sided exact p-value |
| `rank_biserial_synchronized_null_distribution.csv` | Synchronized permutation statistics and evaluable-type counts |
| `Responder_Score_Distribution_AllCells.pdf` | R/NR cosine-score histogram |
| `Reference_Projection_cosine.pdf`, `Reference_Projection_correlation.pdf` | IFNα versus IFNγ projection scatter plots |
| `*_versions.csv` | Package versions recorded for each stage |

Input datasets and generated results are not included in this repository. Local configuration files, inputs, outputs, and intermediate RDS files are excluded through `.gitignore`.
