#!/usr/bin/env Rscript
# Purpose: index original networks or rebuild each cell-type/patient network.
# Inputs: YAML, stage-01 data; archived RDS lists or prepared SCE + scHumanNet.
# Outputs: network_index.rds/csv, network hashes and (build mode only) private RDS.
args_all <- commandArgs(trailingOnly=FALSE)
script <- normalizePath(sub("^--file=", "", args_all[grepl("^--file=",args_all)][1]))
source(file.path(dirname(script),"R","common.R"))

source(file.path(dirname(script),"R","similarity.R"))
source(file.path(dirname(script),"R","construction.R"))
cfg <- cli_config(); record_stage(cfg,"02")
patients <- read_patients(cfg); runs <- read_runs(cfg)
if (cfg$network$mode=="build") {
  load_construction(cfg)
  sce <- readRDS(file.path(cfg$paths$output,"prepared_sce.rds"))
}
index <- list()
for (ri in seq_len(nrow(runs))) {
  run <- runs[ri,,drop=FALSE]; p <- patients_for_run(patients,run)
  if (cfg$network$mode=="build") {
    run_sce <- if (is.na(run$omitted_patient)) sce else sce[,as.character(SummarizedExperiment::colData(sce)[[cfg$metadata$patient]])!=run$omitted_patient]
  }
  for (ct in cfg$cell_type_order) {
    filename <- gsub("{cell_type}",sanitize(ct),cfg$network$filename_template,fixed=TRUE)
    if (cfg$network$mode=="archived") {
      path <- file.path(cfg$paths$archived_networks,filename)
      nets <- readRDS(path)
    } else {
      path <- file.path(cfg$paths$output,"networks",run$run_id,filename)
      if (file.exists(path)) stop("Network already exists; use a new output directory to avoid mixing runs: ",basename(path))
      message(run$run_id," / ",ct)
      sub <- run_sce[,as.character(SummarizedExperiment::colData(run_sce)[[cfg$metadata$fine_cell_type]])==ct]
      labels <- SummarizedExperiment::colData(sub)[[cfg$metadata$patient]]
      if (!all(as.character(labels) %in% p$patient)) stop("Unexpected patient: ",ct)
      nets <- build_networks(sub,labels,cfg)
      dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
      saveRDS(nets,path)
    }
    invisible(edge_vectors(nets,p$patient,cfg$similarity$missing_network_policy))
    index[[length(index)+1L]] <- data.frame(run_id=run$run_id,omitted_patient=run$omitted_patient,cell_type=ct,
      path=normalizePath(path),md5=unname(tools::md5sum(path)))
  }
}
index <- do.call(rbind,index)
saveRDS(index,file.path(cfg$paths$output,"network_index.rds"))
write_csv(index,file.path(cfg$paths$output,"network_index.csv"))
record_stage(cfg,"02")
