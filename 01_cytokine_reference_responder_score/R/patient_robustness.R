# Cell-type-resolved patient-label test on the fixed cosine-score table.
# Source: rank_biserial_global_aggregation_8tests_exact_permutation.R (test 6 only).
# No cell observations or cell-type labels are shuffled independently of patients.
prepare_patient_matrix <- function(dat, cfg) {
  required <- c('donor', 'scvi_fine_CT', 'response', 'Responder_score_cosine')
  if (!all(required %in% names(dat))) stop('Missing projection columns')
  for (key in c('donor', 'scvi_fine_CT')) {
    if (anyNA(dat[[key]]) || any(!nzchar(trimws(dat[[key]])))) stop('Missing patient/cell-type identifiers')
  }
  if (anyDuplicated(paste(dat$donor, dat$scvi_fine_CT, sep='\r'))) stop('Duplicate patient x cell-type rows')
  if (any(!is.finite(dat$Responder_score_cosine))) stop('Non-finite scores')
  dat$group <- ifelse(dat$response == cfg$labels$R, 'R',
                      ifelse(dat$response == cfg$labels$NR, 'NR', NA_character_))
  if (anyNA(dat$group)) stop('Unrecognized or missing response labels')
  groups <- split(dat$group, dat$donor)
  if (any(lengths(lapply(groups, unique)) != 1L)) stop('Conflicting patient response labels')
  patients <- sort(names(groups))
  observed_group <- vapply(patients, function(id) unique(groups[[id]])[[1]], character(1))
  s <- cfg$statistics
  if (length(patients) != s$expected_patients || sum(observed_group=='R') != s$expected_R ||
      sum(observed_group=='NR') != s$expected_NR) stop('Unexpected patient/group counts')
  celltypes <- sort(unique(dat$scvi_fine_CT))
  scores <- matrix(NA_real_, length(patients), length(celltypes), dimnames=list(patients,celltypes))
  scores[cbind(match(dat$donor,patients), match(dat$scvi_fine_CT,celltypes))] <- dat$Responder_score_cosine
  list(scores=scores, group=observed_group)
}

rank_biserial_metrics <- function(scores, ranks, is_R) {
  out <- matrix(NA_real_, ncol(scores), 4L,
                dimnames=list(colnames(scores),c('n_R','n_NR','U_R','r_rb')))
  for (ct in colnames(scores)) {
    valid <- is.finite(scores[,ct])
    r_valid <- valid & is_R
    nr_valid <- valid & !is_R
    n_R <- sum(r_valid); n_NR <- sum(nr_valid)
    out[ct,c('n_R','n_NR')] <- c(n_R,n_NR)
    if (n_R==0L || n_NR==0L) next
    U_R <- sum(ranks[r_valid,ct])-n_R*(n_R+1)/2
    out[ct,c('U_R','r_rb')] <- c(U_R,2*U_R/(n_R*n_NR)-1)
  }
  out
}

aggregate_rank_biserial <- function(metrics, trim_fraction) {
  rrb <- metrics[,'r_rb']
  # Source framework B: omit a cell type only for assignments lacking either group.
  rrb <- rrb[is.finite(rrb)]
  if (!length(rrb)) stop('No evaluable cell types')
  c(n_evaluable=length(rrb), statistic=mean(rrb,trim=trim_fraction))
}

synchronized_rank_biserial_test <- function(prepared, cfg) {
  p <- cfg$rank_biserial
  stopifnot(p$framework=='available_all_celltypes',p$alternative=='greater',
            p$trim_fraction==0.10,p$ties=='average_ranks_half_credit_probability')
  scores <- prepared$scores; groups <- prepared$group
  if (ncol(scores)!=p$expected_celltypes) stop('Unexpected fine-cell-type count')
  ranks <- scores
  for (ct in colnames(scores)) {
    valid <- is.finite(scores[,ct]); values <- scores[valid,ct]
    # The original implementation rejects a non-positive Mann-Whitney variance.
    if (length(unique(values))<2L) stop('Cell type has no score variation')
    ranks[valid,ct] <- rank(values,ties.method='average')
  }
  observed_metrics <- rank_biserial_metrics(scores,ranks,groups=='R')
  observed_aggregate <- aggregate_rank_biserial(observed_metrics,p$trim_fraction)
  if (observed_aggregate['n_evaluable']!=p$expected_observed_evaluable) stop('Unexpected observed evaluability')
  observed <- unname(observed_aggregate['statistic'])
  k <- floor(p$trim_fraction*observed_aggregate['n_evaluable'])
  if (k!=p$expected_trimmed_each_tail) stop('Unexpected observed trim count')
  sorted <- order(observed_metrics[,'r_rb'],na.last=NA)
  retained <- sorted[seq.int(k+1L,length(sorted)-k)]
  trim_status <- rep('not_evaluable',ncol(scores))
  trim_status[sorted] <- 'retained'
  if (k>0) {
    trim_status[head(sorted,k)] <- 'lower_tail'
    trim_status[tail(sorted,k)] <- 'upper_tail'
  }
  assignments <- combn(seq_len(nrow(scores)),cfg$statistics$expected_R)
  n <- ncol(assignments)
  if (n!=cfg$statistics$permutations) stop('Unexpected exact assignment count')
  null <- numeric(n); n_evaluable <- integer(n)
  for (i in seq_len(n)) {
    is_R <- rep(FALSE,nrow(scores)); is_R[assignments[,i]] <- TRUE
    # One shared patient assignment is applied to every fine cell type.
    metrics <- rank_biserial_metrics(scores,ranks,is_R)
    aggregate <- aggregate_rank_biserial(metrics,p$trim_fraction)
    null[i] <- aggregate['statistic']; n_evaluable[i] <- aggregate['n_evaluable']
  }
  if (sum(n_evaluable==23L)!=p$permutations_with_23_celltypes ||
      sum(n_evaluable==22L)!=p$permutations_with_22_celltypes ||
      any(!n_evaluable %in% c(22L,23L))) stop('Unexpected permutation evaluability counts')
  observed_index <- apply(assignments,2L,function(i) identical(i,unname(which(groups=='R'))))
  tolerance <- p$comparison_tolerance*max(1,abs(observed))
  if (sum(observed_index)!=1L || abs(null[observed_index]-observed)>tolerance) stop('Observed assignment check failed')
  extreme <- null>=observed-tolerance
  list(
    celltypes=data.frame(scvi_fine_CT=colnames(scores),patient_coverage=colSums(is.finite(scores)),
      observed_n_R=observed_metrics[,'n_R'],observed_n_NR=observed_metrics[,'n_NR'],
      observed_U_R=observed_metrics[,'U_R'],rank_biserial=observed_metrics[,'r_rb'],
      pairwise_probability_with_half_ties=(observed_metrics[,'r_rb']+1)/2,trim_status=trim_status),
    summary=data.frame(statistic='10% trimmed mean within-cell-type rank-biserial correlation',
      observed_statistic=observed,n_celltypes_observed=observed_aggregate[['n_evaluable']],
      n_trimmed_each_tail=k,n_retained=length(retained),
      pairwise_probability_with_half_ties=(observed+1)/2,
      n_exact_assignments=n,n_as_or_more_extreme=sum(extreme),
      exact_empirical_p_upper_tail=mean(extreme)),
    null=data.frame(permutation_id=seq_len(n),n_evaluable_celltypes=n_evaluable,
      global_statistic=null,as_or_more_extreme=extreme,is_observed_assignment=observed_index)
  )
}
