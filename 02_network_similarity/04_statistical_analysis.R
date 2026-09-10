#!/usr/bin/env Rscript
# Purpose: original NRNR-vs-RR Wilcoxon test, medians and recovered LOPO diagnostics.
# Inputs: stage-03 pairwise_distances.rds.
# Outputs: similarity_summary.csv/rds and optional cell_type_LOPO_robustness_summary.
args_all <- commandArgs(trailingOnly=FALSE)
script <- normalizePath(sub("^--file=", "", args_all[grepl("^--file=",args_all)][1]))
source(file.path(dirname(script),"R","common.R"))

source(file.path(dirname(script),"R","similarity.R"))
source(file.path(dirname(script),"R","robustness.R"))
cfg <- cli_config(); record_stage(cfg,"04")
x <- readRDS(file.path(cfg$paths$output,"pairwise_distances.rds"))
keys <- unique(x[,c("run_id","cell_type")])
s <- do.call(rbind,lapply(seq_len(nrow(keys)),function(i) summarize_pairs(x[x$run_id==keys$run_id[i] & x$cell_type==keys$cell_type[i],])))
rownames(s) <- NULL
saveRDS(s,file.path(cfg$paths$output,"similarity_summary.rds")); write_csv(s,file.path(cfg$paths$output,"similarity_summary.csv"))
if (isTRUE(cfg$robustness$lopo)) {
  r <- summarize_lopo(s,cfg$cell_type_order)
  saveRDS(r,file.path(cfg$paths$output,"cell_type_LOPO_robustness_summary.rds"))
  write_csv(r,file.path(cfg$paths$output,"cell_type_LOPO_robustness_summary.csv"))
}
