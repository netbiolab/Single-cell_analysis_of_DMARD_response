# Call the recovered ACTIONet-HNv3 / SCINET pipeline without redefining algorithms.
# Inputs: prepared SCE, configured reduction/SCINET defaults, HNv3_XC_LLS graph.
# Output: patient-labelled edge lists containing HumanNet LLS and SCINET weights.
load_construction <- function(cfg) {
  suppressPackageStartupMessages({library(ACTIONet);library(SCINET);library(igraph);library(SingleCellExperiment);library(purrr);library(dplyr)})
  if (!is.null(cfg$paths$schumannet_source)) devtools::load_all(cfg$paths$schumannet_source,quiet=TRUE)
  else suppressPackageStartupMessages(library(scHumanNet))
  if (as.character(packageVersion("ACTIONet"))!=cfg$network$expected_actionet || as.character(packageVersion("SCINET"))!=cfg$network$expected_scinet) stop("Construction package version mismatch")
  body_text <- paste(deparse(body(ACTIONet::run.SCINET.clusters)),collapse="\n")
  if (!grepl("HNv3",body_text,fixed=TRUE) || !grepl("log1p",body_text,fixed=TRUE)) stop("Requires modified ACTIONet-HNv3, not stock ACTIONet")
  data("HNv3_XC_LLS",package="scHumanNet",envir=.GlobalEnv)
  if (!exists("graph.hn3",envir=.GlobalEnv)) stop("HNv3_XC_LLS did not supply graph.hn3")
  if (!is.null(cfg$network$global_seed)) set.seed(cfg$network$global_seed)
}
build_networks <- function(sce,labels,cfg) {
  n <- cfg$network
  ace <- ACTIONet::reduce.ace(sce,reduced_dim=n$reduce_dim,max_iter=n$reduce_max_iter,
    assay_name=n$assay,SVD_algorithm=n$svd_algorithm,seed=n$reduce_seed)
  ace[["Labels"]] <- labels
  ace <- ACTIONet::compute.cluster.feature.specificity(ace,ace$Labels,
    "celltype_specificity_scores",assay_name=n$assay)
  graphs <- ACTIONet::run.SCINET.clusters(ace,
    specificity.slot.name="celltype_specificity_scores_feature_specificity",
    min.edge.weight=n$min_edge_weight,spec.sample_no=n$topological_samples,
    thread_no=n$threads,compute.topo.specificity=n$topological_specificity)
  scHumanNet::SortAddLLS(graphs,reference.network=get("graph.hn3",envir=.GlobalEnv))
}
