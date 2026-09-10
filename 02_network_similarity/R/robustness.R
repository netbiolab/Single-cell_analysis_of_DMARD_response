# Recovered LOPO aggregation. Input: full and leave-one-out summary rows.
# Output: per-cell-type median-effect ranges and direction retention.
summarize_lopo <- function(all_similarity, cell_type_order) {
full_similarity <- all_similarity[all_similarity$run_id == "FULL", ]
lopo_similarity <- all_similarity[all_similarity$run_id != "FULL", ]
robustness_rows <- lapply(cell_type_order, function(cell_type) {
  f <- full_similarity[full_similarity$cell_type == cell_type, ]
  x <- lopo_similarity[lopo_similarity$cell_type == cell_type, ]
  same_direction <- sign(x$median_difference_NRNR_minus_RR) ==
    sign(f$median_difference_NRNR_minus_RR)
  data.frame(
    cell_type = cell_type,
    full_median_RR = f$median_RR,
    full_median_R_NR = f$median_R_NR,
    full_median_NRNR = f$median_NRNR,
    full_median_difference_NRNR_minus_RR = f$median_difference_NRNR_minus_RR,
    full_wilcox_p_NRNR_vs_RR = f$wilcox_p_NRNR_vs_RR,
    full_direction = f$direction,
    n_LOPO = nrow(x),
    n_same_median_difference_direction = sum(same_direction),
    proportion_same_median_difference_direction = mean(same_direction),
    LOPO_median_difference_median = median(x$median_difference_NRNR_minus_RR),
    LOPO_median_difference_min = min(x$median_difference_NRNR_minus_RR),
    LOPO_median_difference_max = max(x$median_difference_NRNR_minus_RR),
    max_absolute_LOPO_change_from_full = max(abs(
      x$median_difference_NRNR_minus_RR - f$median_difference_NRNR_minus_RR
    ))
  )
})
robustness_summary <- do.call(rbind, robustness_rows)
rownames(robustness_summary) <- NULL

robustness_summary
}
