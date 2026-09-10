# Network similarity analysis

This directory contains the code used to compare patient-specific functional gene networks between responders (R) and nonresponders (NR) within each cell type.

Networks are inferred using scHumanNet, a HumanNet-compatible ACTIONet implementation, and SCINET. Network similarity is quantified as the Euclidean distance between aligned HumanNet edge-weight vectors, where smaller distances indicate greater network similarity.

The workflow supports:

* analysis of existing patient-specific networks;
* reconstruction of networks from annotated single-cell expression data;
* leave-one-patient-out (LOPO) sensitivity analysis.

## Workflow

```text
Annotated single-cell expression data
        ↓
Cell-type-specific patient grouping
        ↓
ACTIONet feature specificity
        ↓
SCINET network inference on HumanNet v3
        ↓
Edge filtering and HumanNet LLS annotation
        ↓
Patient-network distance calculation
        ↓
RR / R_NR / NRNR comparison
        ↓
Statistical analysis and visualization
```

Existing patient networks can enter directly at the distance-calculation stage after input validation.

---

## Installation

Create the analysis environment using the provided `requirements.txt`:

```bash
conda create -n network_similarity --override-channels \
  -c conda-forge -c bioconda --file requirements.txt

conda activate network_similarity
```

The main network-construction environment uses R 4.0.5, ACTIONet 2.0.18, SCINET 1.0, scHumanNet 0.1.0, SingleCellExperiment 1.12.0, Matrix 1.3-4, and igraph 1.2.6.

For network construction, install the HumanNet-compatible ACTIONet implementation and scHumanNet:

```bash
mkdir -p external

git clone https://github.com/netbiolab/scHumanNet.git external/scHumanNet
git -C external/scHumanNet checkout 79050c292b487abd7ee3d52b27cc57998edced82

R CMD INSTALL external/scHumanNet/packages/ACTIONet_2.0.18_HNv3

# Place a compatible SCINET 1.0 source archive at:
# external/sources/SCINET_1.0.tar.gz

R CMD INSTALL external/sources/SCINET_1.0.tar.gz
R CMD INSTALL external/scHumanNet
```

Preparing Seurat Assay5 input requires a separate R 4.3.1 / Seurat 5.1.0 environment. Analysis of existing patient networks does not require Seurat or network reconstruction.

---

## Inputs and configuration

Copy the example configuration and edit the paths and metadata mapping:

```bash
cp config.example.yaml config.yaml
```

Main inputs are:

| Input              | Description                                                                   |
| ------------------ | ----------------------------------------------------------------------------- |
| Patient manifest   | CSV containing patient identifiers and response groups (`R` or `NR`)          |
| Existing networks  | One RDS file per cell type containing named patient network tables            |
| Expression data    | Annotated Seurat RDS with normalized RNA expression, when rebuilding networks |
| HumanNet resources | HumanNet v3 adjacency and interaction graph supplied through scHumanNet       |

Required metadata fields are configured through:

```yaml
metadata:
  patient: Sample_Name
  fine_cell_type: scvi_fine_CT
  major_cell_type: scvi_major_CT
  response: Response_Final
```

The cell types to analyze and display are specified in `cell_type_order`.

All major network, similarity, statistical, and plotting parameters are defined in `config.example.yaml`.

---

## Network construction

Networks are constructed separately for each patient within each cell type.

For a given cell type:

1. Cells from all patients are processed together.
2. Patient identifiers are used as labels for ACTIONet feature-specificity calculation.
3. `run.SCINET.clusters()` infers patient-specific networks using the HumanNet v3 scaffold.
4. SCINET edges below the configured edge-weight threshold are removed.
5. `scHumanNet::SortAddLLS()` retains HumanNet interactions and adds HumanNet log-likelihood scores (`LLS`) to the retained edges.

A **node** is a gene represented in the HumanNet scaffold and connected by at least one retained interaction.

An **edge** is an undirected HumanNet interaction retained by the SCINET filtering step.

SCINET weights determine whether an edge is retained, whereas the HumanNet `LLS` value is used for the subsequent network-distance calculation.

The pipeline uses the normalized expression supplied by the input object and does not perform additional HVG selection, expression-correlation filtering, or network-size normalization.

---

## Network similarity

For each cell type, the union of all edges observed across available patient networks is constructed.

Each patient's network is represented as an aligned vector of HumanNet LLS values:

```text
v[p,e] = LLS[p,e], if edge e is present in patient p
       = 0,        otherwise
```

The distance between two patient networks is:

```text
D[p,q] = sqrt(sum((v[p,e] - v[q,e])^2))
```

Thus, smaller Euclidean distances indicate more similar network structures.

Three types of patient pairs are evaluated:

* **RR**: responder–responder pairs
* **R_NR**: responder–nonresponder pairs
* **NRNR**: nonresponder–nonresponder pairs

Within-group RR and NRNR pairs retain both pair directions, while each R–NR combination is represented once.

For each cell type, the analysis reports:

* median RR and NRNR distances;
* NRNR-minus-RR median difference;
* NRNR/RR median ratio;
* effect direction;
* two-sided Wilcoxon rank-sum test comparing NRNR and RR distances.

Because pairwise distances share patients and within-group pairs contain both directions, these pairwise observations are not statistically independent patient replicates. The resulting p-values should therefore be interpreted as descriptive comparisons of the pairwise distance distributions.

---

## Missing networks

An absent edge within an available patient network receives an LLS value of zero.

Entirely missing patient networks are controlled by:

```yaml
similarity:
  missing_network_policy: legacy_null_distance_zero
```

`legacy_null_distance_zero` reproduces the original implementation, in which distances involving a missing patient network become zero. These observations are explicitly flagged.

Alternatively:

```yaml
missing_network_policy: error
```

causes the analysis to stop when a required patient network is missing.

A zero generated by the legacy missing-network behavior should not be interpreted as evidence that two networks are identical.

---

## LOPO sensitivity analysis

Leave-one-patient-out analysis evaluates the robustness of the network-similarity pattern.

For each patient:

```text
remove one patient
        ↓
rebuild patient networks
        ↓
recalculate RR / R_NR / NRNR distances
        ↓
recalculate cell-type-level effects
```

The LOPO analysis summarizes:

* consistency of the NRNR-minus-RR effect direction;
* range of effect sizes;
* maximum absolute change relative to the full cohort.

A negative NRNR-minus-RR median difference indicates smaller pairwise network distances among nonresponders.

LOPO is used as a sensitivity analysis; no bootstrap or label-permutation test is performed in this workflow.

---

## Main parameters

Important analysis parameters are defined in `config.example.yaml`.

| Parameter                          |     Default | Description                          |
| ---------------------------------- | ----------: | ------------------------------------ |
| `network.reduce_dim`               |        `50` | ACTIONet reduction dimension         |
| `network.reduce_seed`              |         `0` | ACTIONet reduction seed              |
| `network.min_edge_weight`          |         `2` | Minimum retained SCINET edge weight  |
| `network.topological_samples`      |      `1000` | Topological-specificity samples      |
| `network.threads`                  |         `8` | Network-construction threads         |
| `network.additional_lls_threshold` |      `null` | No additional HumanNet LLS cutoff    |
| `similarity.metric`                | `euclidean` | Distance between aligned LLS vectors |
| `similarity.absent_edge_weight`    |         `0` | Weight assigned to absent edges      |

The original workflow did not set a global random seed. Exact network reconstruction can therefore depend on input order, package implementation, threading, and random-number behavior.

---

## Scripts

| Script                               | Purpose                                                  |
| ------------------------------------ | -------------------------------------------------------- |
| `01_prepare_network_data.R`          | Validate metadata and prepare analysis inputs            |
| `02_build_networks.R`                | Load existing patient networks or construct new networks |
| `03_calculate_similarity.R`          | Calculate pairwise patient-network distances             |
| `04_statistical_analysis.R`          | Summarize distance statistics and LOPO results           |
| `05_plot_results.R`                  | Generate network-distance figures                        |
| `06_build_response_group_networks.R` | Optional pooled response-group network analysis          |
| `R/`                                 | Shared functions                                         |
| `tests/test_numerics.R`              | Synthetic regression tests                               |

---

## Running the main workflow

```bash
cp config.example.yaml config.yaml

# Edit config.yaml before running the analysis.

Rscript 01_prepare_network_data.R --config config.yaml
Rscript 02_build_networks.R --config config.yaml
Rscript 03_calculate_similarity.R --config config.yaml
Rscript 04_statistical_analysis.R --config config.yaml
Rscript 05_plot_results.R --config config.yaml
```

Available configuration profiles are:

| Configuration                  | Purpose                                           |
| ------------------------------ | ------------------------------------------------- |
| `config.example.yaml`          | Analyze existing patient networks                 |
| `configs/rebuild.yaml`         | Reconstruct patient networks from expression data |
| `configs/lopo.yaml`            | Perform full-cohort and LOPO reconstruction       |
| `configs/response_groups.yaml` | Construct optional pooled response-group networks |

Use a separate output directory for each configuration.

The optional response-group network workflow is:

```bash
Rscript 01_prepare_network_data.R --config configs/response_groups.yaml
Rscript 06_build_response_group_networks.R --config configs/response_groups.yaml
```

Synthetic numerical checks can be run with:

```bash
Rscript tests/test_numerics.R
```

---

## Outputs

Stages 03–05 generate the main network-distance results and figures.

Stage 04 additionally generates LOPO robustness summaries when the LOPO configuration is used.

Input datasets, patient-level networks, and generated outputs are not included in this repository. Local configuration files, inputs, outputs, external dependencies, and intermediate networks are excluded through `.gitignore`.
