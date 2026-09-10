#!/usr/bin/env Rscript
# Optional companion: pooled response-by-cell-type networks, separate from similarity.
# Inputs: stage-01 prepared SCE; config response_networks major-cell-type subsets.
# Outputs: per-run pooled networks; original MAIT IRF8 degree-rank summaries.
args_all <- commandArgs(trailingOnly=FALSE)
script <- normalizePath(sub("^--file=", "", args_all[grepl("^--file=",args_all)][1]))
source(file.path(dirname(script),"R","common.R"))


source(file.path(dirname(script),"R","construction.R"))
source(file.path(dirname(script),"R","group_centrality.R"))
cfg <- cli_config(); record_stage(cfg,"06")
if (is.null(cfg$response_networks)) stop("Use configs/response_groups.yaml")
load_construction(cfg)
sce <- readRDS(file.path(cfg$paths$output,"prepared_sce.rds"))
runs <- read_runs(cfg)
results <- list()
for (i in seq_len(nrow(runs))) {
  run <- runs[i,,drop=FALSE]
  run_sce <- if(is.na(run$omitted_patient))sce else sce[,as.character(SummarizedExperiment::colData(sce)[[cfg$metadata$patient]])!=run$omitted_patient]
  for (label in names(cfg$response_networks)) {
    types <- unlist(cfg$response_networks[[label]])
    sub <- run_sce[,as.character(SummarizedExperiment::colData(run_sce)[[cfg$metadata$major_cell_type]]) %in% types]
    if(!ncol(sub)) stop("Empty major-cell-type subset: ",label)
    path <- file.path(cfg$paths$output,"response_group_networks",run$run_id,paste0(label,"_sorted_net_list.rds"))
    if(file.exists(path)) stop("Group network exists; use a new output directory")
    nets <- build_networks(sub,sub$fine_celltype_condition,cfg)
    dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);saveRDS(nets,path)
    if(label=="T_cell") results[[length(results)+1L]] <- irf8_summary(nets,run$run_id,run$omitted_patient,cfg)
  }
}
if(length(results))write_csv(do.call(rbind,results),file.path(cfg$paths$output,"MAIT_IRF8_summary.csv"))
record_stage(cfg,"06")
