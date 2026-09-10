#!/usr/bin/env Rscript
# Input: prepared.rds and frozen centroids or the original cytokine DEG table.
# Output: projection_scores.csv and ordered selected_genes.csv.
# Purpose: calculate the original IFNa-minus-IFNg cosine/Pearson scores.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "R", "common.R"))
cfg <- load_cli()
prepared <- readRDS(out_path(cfg, "prepared.rds"))
if (anyDuplicated(rownames(prepared$expression)) ||
    !identical(colnames(prepared$expression), as.character(prepared$metadata$sample)))
  stop("Prepared expression genes/sample order is invalid")
references <- read_centroids(cfg)
result <- prepared$metadata
for (axis in names(references)) for (method in c("cosine", "correlation")) {
  result[[paste0(axis, "_", method)]] <- apply(prepared$expression, 2, function(x)
    project_sample(x, references[[axis]], cfg$centroid$minimum_shared_genes, method))
}
for (method in c("cosine", "correlation")) {
  result[[paste0("Responder_score_", method)]] <- result[[paste0("IFNa_", method)]] - result[[paste0("IFNg_", method)]]
  result[[paste0("Classification_", method)]] <- ifelse(result[[paste0("Responder_score_", method)]] > 0,
    "Responder (IFNα-like)", "Non-responder (IFNγ-like)")
}
write_csv(result, cfg, "projection_scores.csv")
used <- do.call(rbind, lapply(names(references), function(axis)
  data.frame(axis = axis, rank = seq_along(references[[axis]]), gene = names(references[[axis]]),
    weight = unname(references[[axis]]), present = names(references[[axis]]) %in% rownames(prepared$expression))))
write_csv(used, cfg, "selected_genes.csv")
record_session(cfg, "02")
