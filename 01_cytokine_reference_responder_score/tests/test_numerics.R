#!/usr/bin/env Rscript
# Synthetic tests for projection, aggregation and the synchronized patient-label test.
# No study data or clinical identifiers are included.
script <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
root <- dirname(dirname(normalizePath(script)))
source(file.path(root,'R/common.R'))
source(file.path(root,'R/patient_robustness.R'))
cfg <- yaml::read_yaml(file.path(root,'config.example.yaml'))
assert_error <- function(expr) stopifnot(inherits(tryCatch({force(expr);NULL},error=identity),'error'))
x <- setNames(seq_len(20),paste0('G',seq_len(20)))
for (genes in list(names(x)[1:12],names(x)[9:20])) {
  w <- setNames(seq_along(genes)/2,genes)
  expected <- sum(x[genes]*w)/(sqrt(sum(x[genes]^2))*sqrt(sum(w^2)))
  stopifnot(identical(project_sample(x,w),expected))
}
stopifnot(is.na(project_sample(x,x[1:9])))
raw <- data.frame(cytokine=rep(c('IFN-alpha1','IFN-gamma'),each=5),
 gene=rep(c('A','A','B','C','D'),2),log_fc=rep(c(2,4,1,-1,50),2),
 adj_p_value=rep(c(.01,.01,.01,.01,.05),2))
small_cfg <- cfg;small_cfg$centroid$top_n <- 3
stopifnot(identical(build_centroids(raw,small_cfg)$IFNa,c(A=3,B=1,C=-1)))
expr <- matrix(seq_len(3*21),nrow=3,dimnames=list(c('A','B','C'),paste0('cell',1:21)))
meta <- data.frame(Sample_Name=rep(c('patient_A','patient_B','patient_C'),c(10,10,1)),
 scvi_fine_CT='type_A',Response_Final=rep(c('Good','Bad','Good'),c(10,10,1)),row.names=colnames(expr))
validate_metadata(meta,cfg)
pb <- make_pseudobulk(expr,meta,cfg)
stopifnot(ncol(pb$expression)==2,identical(unname(pb$expression[,1]),unname(rowMeans(expr[,1:10]))))
bad <- meta;bad$Response_Final[1]<-'Bad';assert_error(validate_metadata(bad,cfg))
# All 23 types observed; the sparse type has four patients, two from each group.
scores <- outer(1:15,1:23,function(i,j) ((i*7+j*3) %% 17)+i/10)
rownames(scores) <- sprintf('synthetic_patient_%02d',1:15)
colnames(scores) <- sprintf('type_%02d',1:23)
scores[,23] <- NA_real_;scores[c(1,2,9,10),23] <- c(1,1,2,3)
group <- setNames(c(rep('R',8),rep('NR',7)),rownames(scores))
prepared <- list(scores=scores,group=group)
# Average-rank U must equal the independent sign-of-pairwise-difference definition.
ranks <- apply(scores,2L,rank,na.last='keep',ties.method='average')
metrics <- rank_biserial_metrics(scores,ranks,group=='R')
pairwise <- function(ct,is_R) {
  a <- ct[is_R & is.finite(ct)];b <- ct[!is_R & is.finite(ct)]
  if (!length(a) || !length(b)) return(NA_real_)
  mean(sign(outer(a,b,'-')))
}
expected_effects <- apply(scores,2L,pairwise,is_R=group=='R')
stopifnot(isTRUE(all.equal(unname(metrics[,'r_rb']),unname(expected_effects),tolerance=1e-14)))
result <- synchronized_rank_biserial_test(prepared,cfg)
assignments <- combn(1:15,8)
# Independent all-assignment oracle: recompute pairwise effects and trim after sorting.
expected_null <- apply(assignments,2L,function(ix) {
  effect <- apply(scores,2L,pairwise,is_R=1:15 %in% ix)
  effect <- sort(effect[is.finite(effect)]);k <- floor(length(effect)*.10)
  mean(effect[seq.int(k+1L,length(effect)-k)])
})
stopifnot(isTRUE(all.equal(result$null$global_statistic,unname(expected_null),tolerance=1e-14)),
 sum(result$null$n_evaluable_celltypes==22L)==495L,
 sum(result$null$n_evaluable_celltypes==23L)==5940L,
 result$summary$n_retained==19L,
 result$summary$n_trimmed_each_tail==2L,
 sum(result$null$is_observed_assignment)==1L)
p_expected <- mean(expected_null>=result$summary$observed_statistic-sqrt(.Machine$double.eps))
stopifnot(identical(result$summary$exact_empirical_p_upper_tail,p_expected))
# Projection input validation catches duplicate and conflicting patient records.
dat <- data.frame(donor=rep(rownames(scores),23),scvi_fine_CT=rep(colnames(scores),each=15),
 response=rep(ifelse(group=='R','Good','Bad'),23),Responder_score_cosine=as.vector(scores))
dat <- dat[is.finite(dat$Responder_score_cosine),]
stopifnot(identical(prepare_patient_matrix(dat,cfg)$scores,scores))
assert_error(prepare_patient_matrix(rbind(dat,dat[1,]),cfg))
cat('PASS: centroid selection/projection/means; rank-biserial tie handling; 6435 synchronized nulls; per-assignment eligibility/trimming; one-sided p; input validation\n')
