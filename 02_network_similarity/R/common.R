# Shared I/O and validation. Inputs: YAML and explicitly supplied private files.
# Outputs: validated settings, manifests, and stage-specific run records.
read_config <- function(path) {
  path <- normalizePath(path, mustWork=TRUE)
  cfg <- yaml::read_yaml(path)
  resolve <- function(p) {
    if (is.null(p)) return(NULL)
    if (!grepl("^(/|[A-Za-z]:)", p)) p <- file.path(dirname(path), p)
    normalizePath(p, mustWork=FALSE)
  }
  cfg$paths <- lapply(cfg$paths, resolve)
  fixed <- list(
    c("network", "gene_filter", "intersect_with_HNv3"),
    c("network", "correlation_method", "not_applicable"),
    c("similarity", "metric", "euclidean"),
    c("similarity", "weight", "LLS"),
    c("similarity", "edge_universe", "union_within_cell_type"),
    c("statistics", "test", "wilcox.test"),
    c("statistics", "alternative", "two.sided"),
    c("statistics", "p_adjust", "none"))
  for (x in fixed) if (!identical(cfg[[x[1]]][[x[2]]],x[3])) stop("Unsupported setting: ",paste(x[1:2],collapse="."))
  stopifnot(cfg$network$mode %in% c("archived","build"),
    is.null(cfg$network$additional_lls_threshold),
    isTRUE(cfg$similarity$ordered_within_group_pairs),
    identical(cfg$similarity$normalize,FALSE), cfg$similarity$absent_edge_weight == 0,
    cfg$similarity$missing_network_policy %in% c("legacy_null_distance_zero","error"),
    identical(cfg$statistics$paired,FALSE), is.null(cfg$statistics$exact),
    isTRUE(cfg$statistics$continuity_correction),
    cfg$network$min_edge_weight >= 0, cfg$network$threads >= 1,
    cfg$network$topological_samples >= 1,
    cfg$plot$bubble_size %in% c("legacy_rounded_p","per_cell_type_raw_p"))
  cfg$cell_type_order <- unlist(cfg$cell_type_order, use.names=FALSE)
  if (!length(cfg$cell_type_order) || anyDuplicated(sanitize(cfg$cell_type_order))) stop("Cell types must have unique filename encodings")
  if (cfg$robustness$lopo && cfg$network$mode != "build") stop("LOPO must rebuild networks after patient omission")
  attr(cfg,"config_file") <- path
  cfg
}
cli_config <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  if (length(args)!=2L || args[1]!="--config") stop("Usage: Rscript SCRIPT.R --config CONFIG.yaml")
  read_config(args[2])
}
sanitize <- function(x) gsub("[^A-Za-z0-9]","_",x)
write_csv <- function(x,path) write.csv(x,path,row.names=FALSE,na="")
record_stage <- function(cfg,stage) {
  dir.create(cfg$paths$output,recursive=TRUE,showWarnings=FALSE)
  signature <- cfg; attr(signature,"config_file") <- NULL
  stamp <- file.path(cfg$paths$output,"settings.rds")
  if (file.exists(stamp) && !identical(readRDS(stamp),signature)) stop("Output directory belongs to different settings; choose a new output directory")
  saveRDS(signature,stamp)
  # sessionInfo includes local library paths, so record only package names/versions.
  pkgs <- loadedNamespaces()
  write_csv(data.frame(package=c("R",pkgs), version=c(as.character(getRversion()),vapply(pkgs,function(x)as.character(packageVersion(x)),character(1)))),
            file.path(cfg$paths$output,paste0(stage,"_versions.csv")))
}
read_manifest <- function(cfg) {
  x <- read.csv(cfg$paths$patient_manifest,stringsAsFactors=FALSE,colClasses="character",check.names=FALSE)
  if (!all(c("patient","response_group") %in% names(x))) stop("Manifest needs patient,response_group")
  if (!nrow(x) || anyNA(x$patient) || any(!nzchar(x$patient)) || anyDuplicated(x$patient) || anyNA(x$response_group) || !all(x$response_group %in% c("R","NR"))) stop("Invalid patient manifest")
  if (any(table(factor(x$response_group,levels=c("R","NR")))<2)) stop("At least two patients per group required")
  x[,c("patient","response_group")]
}
read_runs <- function(cfg) readRDS(file.path(cfg$paths$output,"runs.rds"))
read_patients <- function(cfg) readRDS(file.path(cfg$paths$output,"patients.rds"))
patients_for_run <- function(patients,run) {
  if (!is.na(run$omitted_patient)) patients <- patients[patients$patient!=run$omitted_patient,,drop=FALSE]
  if (any(table(factor(patients$response_group,levels=c("R","NR")))<2)) stop("Too few patients after omission")
  patients
}
