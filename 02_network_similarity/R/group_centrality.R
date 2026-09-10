# Recovered MAIT IRF8 diagnostic; not used to compute patient similarity.
# Input: pooled T-cell networks. Output: degree percentile ranks and NR-minus-R difference.
irf8_summary <- function(sorted.net.list,run_id,omitted_patient,cfg) {
  strength.list <- GetCentrality(method = "degree", net.list = sorted.net.list)
  rank.df.final <- CombinePercRank(strength.list)
  r_col <- paste0(cfg$metadata$responder_label,"_MAIT")
  nr_col <- paste0(cfg$metadata$nonresponder_label,"_MAIT")
  if (!all(c(r_col, nr_col) %in% colnames(rank.df.final))) {
    stop("MAIT Good/Bad group networks were not generated for ", run_id)
  }

  diff_pr <- rank.df.final[[nr_col]] - rank.df.final[[r_col]]
  names(diff_pr) <- rownames(rank.df.final)
  if (!("IRF8" %in% names(diff_pr))) stop("IRF8 is absent from the HumanNet gene universe")

  data.frame(
    run_id = run_id,
    omitted_patient = omitted_patient,
    cell_type = "MAIT",
    gene = "IRF8",
    R_degree_percentile_rank = rank.df.final["IRF8", r_col],
    NR_degree_percentile_rank = rank.df.final["IRF8", nr_col],
    differential_centrality_NR_minus_R = diff_pr["IRF8"],
    absolute_differential_hub_rank = rank(-abs(diff_pr), ties.method = "min")["IRF8"],
    NR_greater_than_R_hub_rank = rank(-diff_pr, ties.method = "min")["IRF8"],
    R_degree_hub_rank = rank(
      -setNames(rank.df.final[[r_col]], rownames(rank.df.final)),
      ties.method = "min"
    )["IRF8"],
    NR_degree_hub_rank = rank(
      -setNames(rank.df.final[[nr_col]], rownames(rank.df.final)),
      ties.method = "min"
    )["IRF8"],
    gene_universe_size = length(diff_pr)
  )
}
