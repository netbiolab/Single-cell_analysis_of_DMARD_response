# Environment notes

## Network construction and comparison

Use [requirements.txt](requirements.txt) with the conda command in [README.md](README.md). The specification pins 389 Linux/noarch packages from conda-forge and bioconda, including R 4.0.5 and the dependencies required for network construction, comparison, and plotting.

Three additional R packages must be installed separately:

| Package | Version | Required source |
| --- | --- | --- |
| ACTIONet | 2.0.18 | HumanNet-compatible `packages/ACTIONet_2.0.18_HNv3` from scHumanNet |
| SCINET | 1.0 | Compatible SCINET 1.0 source archive |
| scHumanNet | 0.1.0 | Commit `79050c292b487abd7ee3d52b27cc57998edced82` |

The README provides installation commands. Stock ACTIONet is not equivalent to the HumanNet-compatible implementation. A definitive historical source commit or archive checksum for SCINET was not recovered, so the package version alone does not establish an identical build.

## Seurat input preparation

Preparing Seurat Assay5 input requires a separate R 4.3.1 / Seurat 5.1.0 environment with SingleCellExperiment. Keep this separate from the R 4.0.5 construction environment. Existing patient networks can be analyzed using a patient manifest and network RDS files without Seurat conversion.

## Verification boundary

Package pins were checked against installed conda records and parsed with MatchSpec. Numerical regression checks were performed with the existing construction environment. A clean installation and online availability of all historical builds have not been tested. Compiled source packages also depend on the system/toolchain.
