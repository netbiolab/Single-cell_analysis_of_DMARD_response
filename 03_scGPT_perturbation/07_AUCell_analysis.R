#!/usr/bin/env Rscript
# Purpose: paired SIGLEC1-positive responder AUCell delta signs and random-set null.
# Input: YAML, original Seurat RNA data and decoder prediction H5AD X.
# Output: used gene set, sign fractions, 1000 null values and the combined PDF.
# Source: rerun_manual_ifna_aucell_randomnull_all.R; nofilter/n1000 only.
suppressPackageStartupMessages({
  library(Seurat)
  library(AUCell)
  library(data.table)
  library(ggplot2)
  library(grid)
})

script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(script_file)), "pipeline", "config.R"))
out_dir <- file.path(resolve_path(cfg$paths$output_dir), "AUCell", "nofilter")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
manual_ifna_genes <- unlist(cfg$aucell$gene_set)
filter_label <- "nofilter"
n_iter <- as.integer(cfg$aucell$n_iterations)
if (n_iter != 1000L || cfg$aucell$filter != "nofilter") stop("Only manuscript nofilter/n1000 is supported")
seed_base <- cfg$aucell$seed_base
target_gene <- cfg$aucell$target_gene
if (target_gene != "SIGLEC1") stop("This AUCell figure script is the recorded SIGLEC1 analysis")

trimmed_mean_10 <- function(x) {
  mean(x, trim = cfg$aucell$trim)
}

compose_left_plot <- function(counts_dt, plot_title) {
  path_colors <- c(
    "delta < 0" = "#4C78A8",
    "delta = 0" = "#9D9D9D",
    "delta > 0" = "#F58518"
  )

  counts_dt <- copy(counts_dt)
  counts_dt[, delta_sign := factor(delta_sign, levels = c("delta < 0", "delta = 0", "delta > 0"))]
  setorder(counts_dt, delta_sign)
  counts_dt[, gene_set := factor(gene_set, levels = unique(gene_set))]
  counts_dt[, ypos := cumsum(fraction) - fraction / 2, by = gene_set]

  ggplot(counts_dt, aes(x = gene_set, y = fraction, fill = delta_sign)) +
    geom_col(color = "white", linewidth = 0.3, width = 0.72) +
    geom_text(
      aes(y = ypos, label = ifelse(fraction > 0.025, sprintf("%.1f%%", fraction * 100), "")),
      size = 3.4,
      fontface = "bold"
    ) +
    scale_fill_manual(values = path_colors) +
    scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(round(x * 100), "%")) +
    labs(
      title = plot_title,
      x = "Gene set",
      y = "Fraction of cells",
      fill = "Delta sign"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 15, hjust = 1)
    )
}

compose_right_plot <- function(iter_dt, actual_value, null_mean, p_lower, plot_title) {
  ggplot(iter_dt, aes(x = null_value * 100)) +
    geom_density(fill = "#4C78A8", alpha = 0.30, color = "#2F5D8A", linewidth = 0.8) +
    geom_vline(xintercept = actual_value * 100, color = "firebrick", linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = 0, color = "#7F7F7F", linetype = "dotted", linewidth = 0.6) +
    annotate(
      "label",
      x = Inf,
      y = Inf,
      hjust = 1.02,
      vjust = 1.1,
      size = 3.1,
      label = sprintf(
        "Actual = %.2f%%\nNull mean = %.2f%%\nOne-sided decrease p = %.4f",
        actual_value * 100,
        null_mean * 100,
        p_lower
      )
    ) +
    scale_x_continuous(labels = function(x) sprintf("%.1f%%", x)) +
    labs(
      title = plot_title,
      x = "Trimmed mean delta score (%)",
      y = "Density"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank()
    )
}

save_combined_pdf <- function(p_left, p_right, pdf_file, figure_title) {
  pdf(pdf_file, width = 11, height = 4.8, onefile = TRUE)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2, widths = unit(c(0.9, 1.5), "null"))))
  print(p_left, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_right, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid.text(figure_title, y = unit(0.99, "npc"), gp = gpar(fontsize = 14, fontface = "bold"))
  dev.off()
}

# Direct dense H5AD reader replaces an intermediate Seurat conversion of prediction X.
# This preserves the original converter's pred assay data (which equals prediction X).
library(reticulate)
orig <- readRDS(resolve_path(cfg$paths$original_seurat_rds))
if (DefaultAssay(orig) != "RNA") stop("Historical detection fractions use RNA; original default assay must be RNA")
h5 <- import("h5py", convert = FALSE)
fh <- h5$File(resolve_path(cfg$paths$decoder_predictions), "r")
py_run_string("def read_prediction(f, flag_key):\n import numpy as np\n def strings(ds): return [v.decode() if isinstance(v, bytes) else str(v) for v in ds[()]]\n oi=f['obs'].attrs['_index']; vi=f['var'].attrs['_index']\n flag=f['obs'][flag_key]\n if hasattr(flag, 'keys'):\n  if 'values' not in flag: raise ValueError('Unsupported flag encoding')\n  flags=np.asarray(flag['values'][()], dtype=bool)\n  if 'mask' in flag: flags[np.asarray(flag['mask'][()], dtype=bool)] = False\n else: flags=np.asarray(flag[()], dtype=bool)\n return np.asarray(f['X'][()]), strings(f['obs'][oi]), strings(f['var'][vi]), flags\n")
read_prediction <- py_eval("read_prediction", convert = FALSE)
values <- read_prediction(fh, paste0(target_gene, "_nonzero"))
pred_x <- py_to_r(values[[0L]])
pred_cells <- unlist(py_to_r(values[[1L]]))
pred_genes <- unlist(py_to_r(values[[2L]]))
flag <- as.logical(py_to_r(values[[3L]]))
fh$close()
if (anyDuplicated(pred_cells) || anyDuplicated(pred_genes)) stop("Duplicate prediction axes")
meta <- orig[[]]
m <- cfg$metadata
keep_orig <- which(meta[[m$response]] == m$responder_label & meta[[m$timepoint]] == m$pre_label & meta[[m$cell_type]] == m$ifn_cm_label)
orig_sub <- orig[, keep_orig]
pair_ids <- sub("_[^_]+$", "", pred_cells)
match_idx <- match(pair_ids, colnames(orig_sub))
if (anyNA(match_idx) || anyDuplicated(pair_ids)) stop("Prediction/original pairing failed")
keep <- which(flag)
orig_pair <- orig_sub[, match_idx[keep]]
pert_mat <- t(pred_x[keep, , drop = FALSE])
rownames(pert_mat) <- pred_genes
colnames(pert_mat) <- pred_cells[keep]
common_genes <- intersect(rownames(orig_pair), pred_genes)
orig_mat <- as.matrix(GetAssayData(orig_pair[common_genes, ], assay = DefaultAssay(orig_pair), slot = "data"))
pert_mat <- pert_mat[common_genes, , drop = FALSE]
# Historical code did not seed AUCell's ranking tie breaks. NULL retains that behavior.
if (!is.null(cfg$aucell$ranking_seed)) set.seed(cfg$aucell$ranking_seed)
rankings_orig <- AUCell_buildRankings(orig_mat, plotStats = FALSE, verbose = FALSE)
rankings_pert <- AUCell_buildRankings(pert_mat, plotStats = FALSE, verbose = FALSE)

summary_rows <- list()

pool_genes <- setdiff(common_genes, target_gene)
used_genes <- intersect(manual_ifna_genes, common_genes)
left_title <- "SIGLEC1 AUCell Delta Sign (Nonzero, no filtering)"
right_title_base <- "Random-null density\nAUCell trimmed mean delta, no filtering"

if (!length(used_genes)) stop("No gene-set members in the common gene universe")

  actual_auc_orig <- as.numeric(getAUC(AUCell_calcAUC(list(gs = used_genes), rankings_orig, aucMaxRank = ceiling(nrow(orig_mat) * cfg$aucell$auc_max_rank_fraction), verbose = FALSE))[1, ])
  actual_auc_pert <- as.numeric(getAUC(AUCell_calcAUC(list(gs = used_genes), rankings_pert, aucMaxRank = ceiling(nrow(pert_mat) * cfg$aucell$auc_max_rank_fraction), verbose = FALSE))[1, ])
  actual_delta <- actual_auc_pert - actual_auc_orig
  actual_trimmed_mean <- trimmed_mean_10(actual_delta)

  counts_dt <- data.table(
    gene_set = "manual_IFNa_unique",
    delta_sign = fifelse(actual_delta < 0, "delta < 0", fifelse(actual_delta > 0, "delta > 0", "delta = 0"))
  )[, .N, by = .(gene_set, delta_sign)]
  counts_dt[, fraction := N / sum(N)]

  fwrite(counts_dt, file.path(out_dir, paste0("delta_sign_", filter_label, ".csv")))
  fwrite(data.table(gene = used_genes), file.path(out_dir, paste0("used_genes_", filter_label, ".csv")))
  left_plot <- compose_left_plot(counts_dt[, .(gene_set, delta_sign, fraction)], left_title)

    set.seed(seed_base)
    null_vals <- numeric(n_iter)
    for (i in seq_len(n_iter)) {
      rand_genes <- sample(pool_genes, length(used_genes), replace = FALSE)
      auc_orig_rand <- as.numeric(getAUC(AUCell_calcAUC(list(gs = rand_genes), rankings_orig, aucMaxRank = ceiling(nrow(orig_mat) * cfg$aucell$auc_max_rank_fraction), verbose = FALSE))[1, ])
      auc_pert_rand <- as.numeric(getAUC(AUCell_calcAUC(list(gs = rand_genes), rankings_pert, aucMaxRank = ceiling(nrow(pert_mat) * cfg$aucell$auc_max_rank_fraction), verbose = FALSE))[1, ])
      null_vals[i] <- trimmed_mean_10(auc_pert_rand - auc_orig_rand)
    }

    p_lower <- mean(null_vals <= actual_trimmed_mean)
    null_mean <- mean(null_vals)
    iter_dt <- data.table(
      metric = "AUCell",
      gene_set = "manual_IFNa_unique",
      filter_label = filter_label,
      feature = "trimmed_mean",
      iteration = seq_len(n_iter),
      null_value = null_vals
    )

    right_plot <- compose_right_plot(iter_dt, actual_trimmed_mean, null_mean, p_lower, right_title_base)

    prefix <- sprintf("manual_IFNa_unique_AUCell_%s_trimmed_mean", filter_label)
    pdf_file <- file.path(out_dir, sprintf("%s_combined_siglec1_nonzero_n%s.pdf", prefix, n_iter))
    csv_file <- file.path(out_dir, sprintf("%s_null_iterations_n%s.csv", prefix, n_iter))

    fwrite(iter_dt, csv_file)
    save_combined_pdf(
      left_plot,
      right_plot,
      pdf_file,
      "manual_IFNa_unique AUCell Summary, SIGLEC1 nonzero"
    )

    summary_rows[[paste(filter_label, n_iter, sep = "_")]] <- data.table(
      metric = "AUCell",
      gene_set = "manual_IFNa_unique",
      filter_label = filter_label,
      feature = "trimmed_mean",
      n_used_genes = length(used_genes),
      n_pool_genes = length(pool_genes),
      n_iter = n_iter,
      actual_value = actual_trimmed_mean,
      null_mean = null_mean,
      null_sd = sd(null_vals),
      empirical_p_lower = p_lower,
      empirical_p_upper = mean(null_vals >= actual_trimmed_mean),
      two_sided_empirical_p = min(1, 2 * min(mean(null_vals <= actual_trimmed_mean), mean(null_vals >= actual_trimmed_mean))),
      output_pdf = pdf_file,
      output_csv = csv_file
    )

summary_dt <- rbindlist(summary_rows)
fwrite(summary_dt, file.path(out_dir, "manual_IFNa_unique_AUCell_nofilter_trimmed_mean_n1000_summary.csv"))
cat("Saved outputs to:\n", out_dir, "\n", sep = "")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
