#!/usr/bin/env Rscript
# Input: preprocessed Seurat RDS with RNA/data and configured metadata.
# Output: prepared.rds (donor x cell-type means).
# Purpose: preserve the source expression scale, cell order and aggregation scope.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "R", "common.R"))
cfg <- load_cli()
object <- readRDS(cfg$paths$seurat)
validate_metadata(object@meta.data, cfg)
expr <- Seurat::GetAssayData(object, assay = "RNA", slot = "data")
prepared <- make_pseudobulk(expr, object@meta.data, cfg)
saveRDS(prepared, out_path(cfg, "prepared.rds"))
record_session(cfg, "01")
