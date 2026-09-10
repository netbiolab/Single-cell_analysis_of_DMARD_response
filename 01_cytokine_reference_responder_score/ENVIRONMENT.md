# Environment notes

Use the conda installation command in [README.md](README.md). The Linux requirements pin 255 public conda-forge packages by version and build, including the numerical, rendering, and runtime dependencies of the analysis.

The target uses R 4.3.1, Seurat 5.1.0, SeuratObject 5.0.2, Matrix 1.6-5, dplyr 1.1.4, ggplot2 3.5.1, and yaml 2.3.7. Conda represents Matrix 1.6-5 as `1.6_5`; R may display `1.6.5`.

Matrix 1.6-5 replaces an earlier R-Forge development build labelled 1.7-0. This compatibility change was checked without changing the scoring or statistical implementation.

A fresh environment was created from cached public package archives, and all installed version/build records matched the requirements. Required namespaces and Cairo PDF support loaded successfully. The four-stage workflow and synthetic tests passed. Prepared expression profiles, scores, statistical summaries, and complete exact-permutation distributions matched the earlier validated outputs. PDF generation succeeded; byte-identical rendering is not claimed. Public-channel download availability was not separately tested.
