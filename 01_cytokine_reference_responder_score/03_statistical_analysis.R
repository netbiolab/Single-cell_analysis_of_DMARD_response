#!/usr/bin/env Rscript
# Input: projection_scores.csv from stage 02.
# Output: pooled KS/group summaries, patient-median exact test, and the manuscript
# 10% trimmed rank-biserial synchronized patient-label permutation test.
script <- sub('^--file=', '', grep('^--file=', commandArgs(), value=TRUE)[1])
script_dir <- dirname(normalizePath(script))
source(file.path(script_dir,'R','common.R'))
source(file.path(script_dir,'R','patient_robustness.R'))
cfg <- load_cli()
dat <- read.csv(out_path(cfg,'projection_scores.csv'),check.names=FALSE)
prepared <- prepare_patient_matrix(dat,cfg)
if (isTRUE(cfg$statistics$group_comparison)) {
  r <- dat$Responder_score_cosine[dat$response==cfg$labels$R]
  nr <- dat$Responder_score_cosine[dat$response==cfg$labels$NR]
  test <- ks.test(r,nr,alternative='two.sided',exact=NULL)
  write_csv(data.frame(test=test$method,unit='donor_x_celltype',
    statistic=unname(test$statistic),p_value=test$p.value,adjustment='none'),cfg,'statistical_tests.csv')
  write_csv(do.call(rbind,lapply(c(cfg$labels$R,cfg$labels$NR),function(group) {
    x <- dat$Responder_score_cosine[dat$response==group]
    data.frame(response=group,n=length(x),mean=mean(x),median=median(x),sd=sd(x))
  })),cfg,'group_summary.csv')
}
if (isTRUE(cfg$statistics$patient_exact_permutation)) {
  result <- patient_median_test(prepared,cfg)
  write_csv(result$patients,cfg,'robustness_patient_level_scores.csv')
  write_csv(result$summary,cfg,'robustness_exact_permutation_summary.csv')
  write_csv(result$null,cfg,'robustness_exact_permutation_null_distribution.csv')
}
if (isTRUE(cfg$statistics$rank_biserial_permutation)) {
  result <- synchronized_rank_biserial_test(prepared,cfg)
  write_csv(result$celltypes,cfg,'rank_biserial_celltype_effects.csv')
  write_csv(result$summary,cfg,'rank_biserial_trimmed_mean_exact_summary.csv')
  write_csv(result$null,cfg,'rank_biserial_synchronized_null_distribution.csv')
}
record_session(cfg,'03')
