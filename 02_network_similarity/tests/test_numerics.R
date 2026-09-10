#!/usr/bin/env Rscript
# Synthetic regression against the original literal edge-union and nested loops.
# No study data are needed. Run from any working directory.
a <- commandArgs(FALSE); here <- dirname(normalizePath(sub("^--file=","",a[grepl("^--file=",a)][1])))
source(file.path(here,"..","R","similarity.R"))
source(file.path(here,"..","R","common.R"))
net <- function(a,b,w) data.frame(gene1=a,gene2=b,LLS=w,scinet_weight=100-w)
nets <- list(R1=net(c("B","C"),c("A","B"),c(2,3)),R2=net("A","B",2),
             N1=net("D","A",4),N2=net(c("B","D"),c("C","A"),c(3,4)))
patients <- data.frame(patient=c("R1","R2","N1","N2"),response_group=c("R","R","NR","NR"))
# Literal source implementation (independent, deliberately slow).
original <- function(sorted.net.list,good,bad) {
 all.edge <- c()
 for(i in 1:length(sorted.net.list)) {
  net <- sorted.net.list[[i]]
  net.edge <- unname(unlist(apply(net,1,function(x)paste0(sort(c(x[1],x[2])),collapse="|"))))
  net[["edge"]] <- net.edge; sorted.net.list[[i]] <- net; all.edge <- c(all.edge,net.edge)
 }
 all.edge <- unique(all.edge); distance.list <- list()
 for(i in 1:length(sorted.net.list)) {
  net <- sorted.net.list[[i]]; vector <- c()
  for(j in 1:length(all.edge)) {
   if(all.edge[j] %in% net$edge) vector <- c(vector,net[net$edge==all.edge[j],]$LLS)
   else vector <- c(vector,0)
  }
  distance.list[[i]] <- vector
 }
 names(distance.list) <- names(sorted.net.list)
 euclidean <- function(a,b)sqrt(sum((a-b)^2))
 NRR <- c();RR <- c();NRNR <- c()
 for(i in 1:length(good))for(j in 1:length(bad))NRR<-c(NRR,euclidean(distance.list[[good[i]]],distance.list[[bad[j]]]))
 for(i in 1:length(good))for(j in 1:length(good))if(i!=j)RR<-c(RR,euclidean(distance.list[[good[i]]],distance.list[[good[j]]]))
 for(i in 1:length(bad))for(j in 1:length(bad))if(i!=j)NRNR<-c(NRNR,euclidean(distance.list[[bad[i]]],distance.list[[bad[j]]]))
 c(RR,NRR,NRNR)
}
x <- calculate_pairs(nets,patients)
stopifnot(identical(x$distance,original(nets,c("R1","R2"),c("N1","N2"))),
          identical(as.integer(table(x$group)),c(2L,4L,2L)),
          all(x$distance[x$group=="RR"]==3),!any(x$missing_patient_network))
# Missing whole network must retain the source's NULL subtraction -> zero.
y <- calculate_pairs(nets[-4],patients)
stopifnot(identical(y$distance,original(nets[-4],c("R1","R2"),c("N1","N2"))),
          all(y$distance[y$missing_patient_network]==0))
fails <- function(expr) inherits(tryCatch({force(expr);NULL},error=function(e)e),"error")
stopifnot(fails(calculate_pairs(nets[-4],patients,missing_network_policy="error")))
bad <- nets;bad[[1]]<-rbind(bad[[1]],bad[[1]][1,]);stopifnot(fails(edge_vectors(bad,patients$patient)))
bad <- nets;bad[[1]]$LLS[1]<-NA;stopifnot(fails(edge_vectors(bad,patients$patient)))
z <- summarize_pairs(x)
stopifnot(identical(z$wilcox_p_NRNR_vs_RR,suppressWarnings(wilcox.test(x$distance[x$group=="NRNR"],x$distance[x$group=="RR"])$p.value)))
cfg <- read_config(file.path(here,"..","config.example.yaml"));stopifnot(cfg$network$min_edge_weight==2)
cat("PASS: literal-source parity, edge reversal/union/zero-fill, pair multiplicity, missing-network behavior, invalid input, Wilcoxon.\n")
# Requested release scope: retained statistical schema and LOPO median diagnostics.
stopifnot(identical(names(z),c("run_id","omitted_patient","cell_type","median_RR",
  "median_R_NR","median_NRNR","median_difference_NRNR_minus_RR",
  "median_ratio_NRNR_over_RR","wilcox_p_NRNR_vs_RR","direction")))
source(file.path(here,"..","R","robustness.R"))
full <- z; full$median_difference_NRNR_minus_RR <- 2
leave_a <- full; leave_a$run_id <- "LOPO_A"; leave_a$median_difference_NRNR_minus_RR <- 1
leave_b <- full; leave_b$run_id <- "LOPO_B"; leave_b$median_difference_NRNR_minus_RR <- -1
r <- summarize_lopo(rbind(full,leave_a,leave_b),full$cell_type)
stopifnot(!any(grepl("cliffs|most_influential",names(r))),
  r$n_LOPO==2,r$proportion_same_median_difference_direction==0.5,
  r$LOPO_median_difference_median==0,r$max_absolute_LOPO_change_from_full==3)
cat("PASS: reduced summary schema and retained LOPO median diagnostics.\n")
