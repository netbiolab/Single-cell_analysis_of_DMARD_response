#!/usr/bin/env Rscript
# Purpose: validate revised patient groups; optionally convert annotated Seurat to SCE.
# Inputs: YAML, patient_manifest.csv, optional Seurat RDS with existing RNA/data.
# Outputs: patients.rds, runs.rds, optional prepared_sce.rds and cell_counts.csv.
args_all <- commandArgs(trailingOnly=FALSE)
script <- normalizePath(sub("^--file=", "", args_all[grepl("^--file=",args_all)][1]))
source(file.path(dirname(script),"R","common.R"))

cfg <- cli_config()
patients <- read_manifest(cfg)
record_stage(cfg,"01")
if (!is.null(cfg$paths$seurat)) {
  suppressPackageStartupMessages({library(Seurat);library(SingleCellExperiment)})
  obj <- readRDS(cfg$paths$seurat)
  cols <- unlist(cfg$metadata[c("patient","fine_cell_type","major_cell_type")])
  if (!all(cols %in% names(obj[[]]))) stop("Missing required Seurat metadata columns")
  md <- obj[[]]
  if (anyNA(md[,cols,drop=FALSE])) stop("Missing patient/cell-type metadata")
  ids <- as.character(md[[cfg$metadata$patient]])
  if (!setequal(unique(ids),patients$patient)) stop("Seurat and manifest patient sets differ")
  if (!all(cfg$cell_type_order %in% as.character(md[[cfg$metadata$fine_cell_type]]))) stop("Missing requested cell types")
  # The clinical revision is explicit private input, replacing original hardcoded IDs.
  groups <- patients$response_group[match(ids,patients$patient)]
  obj[[cfg$metadata$response]] <- ifelse(groups=="R",cfg$metadata$responder_label,cfg$metadata$nonresponder_label)
  obj$fine_celltype_condition <- paste(obj[[]][[cfg$metadata$response]],md[[cfg$metadata$fine_cell_type]],sep="_")
  obj$major_celltype_condition <- paste(obj[[]][[cfg$metadata$response]],md[[cfg$metadata$major_cell_type]],sep="_")
  obj[["RNA"]] <- as(obj[["RNA"]],Class="Assay")
  # Match the original conversion; fail if its active assay would select another scale.
  if (DefaultAssay(obj)!="RNA") stop("Original workflow requires active RNA assay")
  sce <- Seurat::as.SingleCellExperiment(obj)
  if (!"logcounts" %in% SummarizedExperiment::assayNames(sce)) stop("No normalized logcounts assay")
  saveRDS(sce,file.path(cfg$paths$output,"prepared_sce.rds"))
  counts <- as.data.frame(table(patient=ids,cell_type=md[[cfg$metadata$fine_cell_type]]))
  names(counts)[3] <- "n_cells"
  write_csv(counts,file.path(cfg$paths$output,"cell_counts.csv"))
} else if (cfg$network$mode=="build") stop("Building requires paths.seurat")
# Keep manifest order; original pair loops use responder order followed by NR order.
patients <- rbind(patients[patients$response_group=="R",],patients[patients$response_group=="NR",])
runs <- data.frame(run_id="FULL",omitted_patient=NA_character_,stringsAsFactors=FALSE)
if (isTRUE(cfg$robustness$lopo)) runs <- rbind(runs,data.frame(run_id=paste0("LOPO_",sanitize(patients$patient)),omitted_patient=patients$patient))
if (anyDuplicated(runs$run_id)) stop("Patient IDs collide after filename sanitization")
saveRDS(patients,file.path(cfg$paths$output,"patients.rds"))
saveRDS(runs,file.path(cfg$paths$output,"runs.rds"))
write_csv(runs,file.path(cfg$paths$output,"run_manifest.csv"))
record_stage(cfg,"01")
