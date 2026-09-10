#!/usr/bin/env Rscript
# Purpose: align the edge union, zero-fill missing edges, calculate Euclidean distances.
# Inputs: stage-02 network index and patient manifest; verified network RDS lists.
# Outputs: pairwise_distances.rds/csv.gz, with original ordered RR/NRNR pairs.
args_all <- commandArgs(trailingOnly=FALSE)
script <- normalizePath(sub("^--file=", "", args_all[grepl("^--file=",args_all)][1]))
source(file.path(dirname(script),"R","common.R"))

source(file.path(dirname(script),"R","similarity.R"))
cfg <- cli_config(); record_stage(cfg,"03")
index <- readRDS(file.path(cfg$paths$output,"network_index.rds")); patients <- read_patients(cfg)
results <- lapply(seq_len(nrow(index)),function(i) {
  x <- index[i,,drop=FALSE]
  if (unname(tools::md5sum(x$path))!=x$md5) stop("Network changed after indexing")
  p <- patients_for_run(patients,x)
  calculate_pairs(readRDS(x$path),p,x$run_id,x$omitted_patient,x$cell_type,cfg$similarity$missing_network_policy)
})
x <- do.call(rbind,results); rownames(x) <- NULL
saveRDS(x,file.path(cfg$paths$output,"pairwise_distances.rds"))
con <- gzfile(file.path(cfg$paths$output,"pairwise_distances.csv.gz"),"wt")
write_csv(x,con); close(con)
