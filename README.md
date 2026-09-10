# RA-DMARD-response

## Overview

Analysis code for studying disease-modifying antirheumatic drug (DMARD) treatment response in rheumatoid arthritis using single-cell RNA sequencing of pretreatment peripheral blood mononuclear cells (PBMCs). The analyses compare responders (R) and nonresponders (NR).

## Main analyses

- **[Responder Score](01_cytokine_reference_responder_score/README.md):** cytokine-reference-based scoring and comparison of response groups.
- **[Network similarity](02_network_similarity/README.md):** comparison of patient-specific functional gene networks within cell types.
- **[scGPT in silico perturbation](03_scGPT_perturbation/README.md):** simulated gene-expression perturbations and evaluation of their effects in IFN CM cells.

Each analysis README describes its methods, inputs, and execution steps.

## Repository structure

```text
RA-DMARD-response/
├── README.md
├── .gitignore
├── 01_cytokine_reference_responder_score/
├── 02_network_similarity/
└── 03_scGPT_perturbation/
```

## Data availability

Raw and processed single-cell data are provided through the Gene Expression Omnibus (GEO) and are not included in this repository. Patient-identifiable information and private clinical data are not distributed here.

## Installation

Environments are analysis-specific. Each analysis directory provides a `requirements.txt` for conda. The scGPT directory additionally provides `requirements-pip.txt` and `requirements-r.txt` for its Python dependencies and separate R environment.

Follow the installation instructions in the corresponding analysis README, including the additional scHumanNet/ACTIONet/SCINET and scGPT/flash-attention setup where applicable. Network construction and Seurat input preparation also use separate environments.

## Usage

Clone the repository and enter the analysis directory of interest. For example, to run the Responder Score workflow, replace `<repository-url>` with this repository's clone URL:

```bash
git clone <repository-url> RA-DMARD-response
cd RA-DMARD-response/01_cytokine_reference_responder_score
cp config.example.yaml config.yaml
```

Read the analysis README, prepare its environment and required inputs, then edit `config.yaml` and follow the documented execution order.

## Reproducibility

Each analysis records major parameters and input/output paths in `config.example.yaml`; package requirements and numerical tests are provided alongside the code. Large intermediate files, model checkpoints, and generated outputs are excluded from GitHub through `.gitignore`.
