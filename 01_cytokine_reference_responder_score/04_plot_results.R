#!/usr/bin/env Rscript
# Input: stage-02 and stage-03 CSVs. Output: distribution and projection PDFs.
# Purpose: render recovered result types without recalculating or changing scores.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(normalizePath(script)), "R", "common.R"))
suppressPackageStartupMessages(library(ggplot2))
cfg <- load_cli()
colors <- setNames(unlist(cfg$plot$colors), c(cfg$labels$R, cfg$labels$NR))
dat <- read.csv(out_path(cfg, "projection_scores.csv"))
test <- read.csv(out_path(cfg, "statistical_tests.csv"))
p <- ggplot(dat, aes(Responder_score_cosine, fill = response)) +
  geom_histogram(bins = cfg$plot$bins, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed") + scale_fill_manual(values = colors) +
  labs(x = "Responder Score (IFNα - IFNγ)", y = "Number of Samples",
    subtitle = paste0("Donor × cell type | KS p = ", signif(test$p_value, 3))) + theme_bw()
ggsave(device = grDevices::cairo_pdf, filename = out_path(cfg, "Responder_Score_Distribution_AllCells.pdf"), p, width = 10, height = 6)
for (method in c("cosine", "correlation")) {
  p <- ggplot(dat, aes(.data[[paste0("IFNa_", method)]], .data[[paste0("IFNg_", method)]], color = response)) +
    geom_point(size = 3, alpha = 0.7) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = colors) + labs(x = paste("IFNα", method), y = paste("IFNγ", method)) + theme_bw()
  ggsave(device = grDevices::cairo_pdf, filename = out_path(cfg, paste0("Reference_Projection_", method, ".pdf")), p, width = 10, height = 8)
}
record_session(cfg, "04")
