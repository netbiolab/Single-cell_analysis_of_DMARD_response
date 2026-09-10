# Edge alignment and original ordered patient-pair Euclidean distances.
# Input: named list of edge tables (endpoints in first 2 columns; numeric LLS).
# Output: long pair table. Invalid/duplicate edges fail rather than recycle vectors.
edge_vectors <- function(nets, patients, missing_network_policy="legacy_null_distance_zero") {
  if (!is.list(nets) || is.null(names(nets)) || anyDuplicated(names(nets)) || !all(names(nets) %in% patients)) stop("Network names must be a subset of the supplied patient set")
  if (missing_network_policy=="error" && !setequal(names(nets),patients)) stop("Missing patient network")
  edges <- lapply(nets,function(n) {
    if (!is.data.frame(n) || nrow(n)==0L || ncol(n)<3L || !"LLS" %in% names(n) || !is.numeric(n$LLS) || any(!is.finite(n$LLS))) stop("Network requires nonempty edge table and finite numeric LLS")
    endpoints <- n[,1:2,drop=FALSE]
    if (anyNA(endpoints) || any(!nzchar(as.matrix(endpoints))) || any(grepl("|",as.matrix(endpoints),fixed=TRUE))) stop("Invalid endpoint identifier")
    e <- unname(apply(endpoints,1,function(x)paste0(sort(as.character(x)),collapse="|")))
    if (anyDuplicated(e)) stop("Duplicate undirected edge; source would misalign vectors")
    e
  })
  universe <- unique(unlist(edges,use.names=FALSE))
  vectors <- lapply(seq_along(nets),function(i) {
    idx <- match(universe,edges[[i]])
    v <- numeric(length(universe)); found <- !is.na(idx)
    v[found] <- nets[[i]]$LLS[idx[found]]
    v
  })
  names(vectors) <- names(nets)
  vectors
}
calculate_pairs <- function(nets,patients,run_id="FULL",omitted_patient=NA_character_,cell_type="synthetic",missing_network_policy="legacy_null_distance_zero") {
  v <- edge_vectors(nets,patients$patient,missing_network_policy)
  good <- patients$patient[patients$response_group=="R"]
  bad <- patients$patient[patients$response_group=="NR"]
  pairs <- function(a,b,group,within=FALSE) {
    do.call(rbind,lapply(a,function(i) {
      js <- if(within)b[b!=i] else b
      data.frame(group=group,patient_1=i,patient_2=js,
        distance=vapply(js,function(j)sqrt(sum((v[[i]]-v[[j]])^2)),numeric(1)),stringsAsFactors=FALSE)
    }))
  }
  x <- rbind(pairs(good,good,"RR",TRUE),pairs(good,bad,"R_NR"),pairs(bad,bad,"NRNR",TRUE))
  x$missing_patient_network <- !x$patient_1 %in% names(nets) | !x$patient_2 %in% names(nets)
  x$run_id <- run_id; x$omitted_patient <- omitted_patient; x$cell_type <- cell_type
  rownames(x) <- NULL
  x[,c("run_id","omitted_patient","cell_type","group","patient_1","patient_2","distance","missing_patient_network")]
}
summarize_pairs <- function(x) {
  rr <- x$distance[x$group=="RR"]; rn <- x$distance[x$group=="R_NR"]; nn <- x$distance[x$group=="NRNR"]
  if (!length(rr)||!length(nn)||!length(rn)||any(!is.finite(x$distance))) stop("Incomplete pair table")
  difference <- median(nn)-median(rr)
  data.frame(run_id=x$run_id[1],omitted_patient=x$omitted_patient[1],cell_type=x$cell_type[1],
    median_RR=median(rr),median_R_NR=median(rn),median_NRNR=median(nn),
    median_difference_NRNR_minus_RR=difference,median_ratio_NRNR_over_RR=median(nn)/median(rr),
    wilcox_p_NRNR_vs_RR=suppressWarnings(wilcox.test(nn,rr,alternative="two.sided",paired=FALSE,exact=NULL,correct=TRUE)$p.value),
    direction=if(difference<0)"NR_convergence" else if(difference>0)"R_convergence" else "no_difference")
}
