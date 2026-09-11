# Shared configuration, gene selection and numerical routines.
# Source from the numbered entry points; paths resolve relative to the YAML file.
read_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  root <- dirname(normalizePath(path, mustWork = TRUE))
  for (key in names(cfg$paths)) {
    value <- cfg$paths[[key]]
    if (!is.null(value) && !grepl("^(/|[A-Za-z]:)", value)) value <- file.path(root, value)
    cfg$paths[[key]] <- value
  }
  stopifnot(cfg$workflow == "centroid")
  stopifnot(cfg$centroid$expression_aggregation == "mean_RNA_data",
    cfg$centroid$weights == "mean_log_fc", !cfg$centroid$positive_logfc_filter,
    !cfg$centroid$zscore, cfg$statistics$p_adjust == "none")
  dir.create(cfg$paths$output, recursive = TRUE, showWarnings = FALSE)
  cfg
}

load_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L || args[1] != "--config") stop("Usage: Rscript SCRIPT.R --config CONFIG.yaml")
  read_config(args[2])
}

out_path <- function(cfg, name) file.path(cfg$paths$output, name)
write_csv <- function(x, cfg, name) write.csv(x, out_path(cfg, name), row.names = FALSE)
record_session <- function(cfg, stage) {
  # Deliberately omit local library paths and host/user names from this record.
  packages <- c("Seurat", "SeuratObject", "Matrix", "dplyr", "yaml", "ggplot2")
  versions <- vapply(packages, function(x)
    tryCatch(as.character(utils::packageVersion(x)), error = function(e) NA_character_), character(1))
  write_csv(data.frame(package = c("R", packages), version = c(as.character(getRversion()), versions)),
            cfg, paste0(stage, "_versions.csv"))
}

validate_metadata <- function(meta, cfg) {
  cols <- cfg$columns
  needed <- c(cols$patient, cols$celltype, cols$response)
  if (!all(needed %in% names(meta))) stop("Missing required metadata columns: ", paste(setdiff(needed, names(meta)), collapse = ", "))
  for (key in c(cols$patient, cols$celltype)) {
    if (anyNA(meta[[key]]) || any(!nzchar(as.character(meta[[key]])))) stop("Missing patient/cell-type metadata")
  }
  labels <- as.character(meta[[cols$response]])
  if (any(!is.na(labels) & !labels %in% c(cfg$labels$R, cfg$labels$NR))) stop("Unexpected response labels")
  groups <- split(labels, meta[[cols$patient]])
  if (any(vapply(groups, function(x) length(unique(x[!is.na(x)])) != 1L, logical(1)))) stop("Patient response labels are missing or inconsistent")
}

build_centroids <- function(deg, cfg) {
  c <- cfg$centroid
  if (!all(c("cytokine", "adj_p_value", "gene", "log_fc") %in% names(deg))) stop("Invalid DEG table")
  result <- lapply(c("IFNa", "IFNg"), function(axis) {
    ranked <- deg |>
      dplyr::filter(.data$cytokine %in% unlist(c[[axis]]), .data$adj_p_value < c$adjusted_p_cutoff) |>
      dplyr::group_by(.data$gene) |>
      dplyr::summarize(mean_log_fc = mean(.data$log_fc), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(.data$mean_log_fc))
    ranked <- head(ranked, c$top_n)
    setNames(ranked$mean_log_fc, ranked$gene)
  })
  setNames(result, c("IFNa", "IFNg"))
}

read_centroids <- function(cfg) {
  if (cfg$centroid$source == "deg_table") return(build_centroids(read.csv(cfg$paths$cytokine_degs), cfg))
  stopifnot(cfg$centroid$source == "frozen")
  tab <- read.delim(cfg$paths$centroids, check.names = FALSE)
  setNames(lapply(c("IFNa", "IFNg"), function(axis) {
    x <- tab[tab$axis == axis, ]; x <- x[order(x$rank), ]
    if (nrow(x) < cfg$centroid$top_n) stop("Frozen reference is too short for requested top_n")
    x <- head(x, cfg$centroid$top_n)
    setNames(x$weight, x$gene)
  }), c("IFNa", "IFNg"))
}

project_sample <- function(expr, centroid, minimum = 10L, method = "cosine") {
  genes <- intersect(names(expr), names(centroid))
  if (length(genes) < minimum) return(NA_real_)
  x <- expr[genes]; y <- centroid[genes]
  if (method == "cosine") sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
  else cor(x, y, method = "pearson")
}

make_pseudobulk <- function(expr, meta, cfg) {
  cols <- cfg$columns
  donors <- unique(meta[[cols$patient]])
  celltypes <- unique(meta[[cols$celltype]])
  profiles <- list(); rows <- list()
  for (donor in donors) for (ct in celltypes) {
    cells <- rownames(meta)[meta[[cols$patient]] == donor & meta[[cols$celltype]] == ct]
    if (length(cells) < cfg$centroid$minimum_cells) next
    id <- paste0(donor, "_", ct)
    if (id %in% names(profiles)) stop("Ambiguous sample IDs")
    x <- expr[, cells, drop = FALSE]
    profiles[[id]] <- if (inherits(x, "sparseMatrix")) Matrix::rowMeans(x) else rowMeans(x)
    rows[[id]] <- data.frame(sample = id, donor = as.character(donor), scvi_fine_CT = as.character(ct),
      n_cells = as.character(length(cells)), response = as.character(unique(meta[cells, cols$response])[1]))
  }
  if (!length(profiles)) stop("No groups pass minimum_cells")
  list(expression = do.call(cbind, profiles), metadata = do.call(rbind, rows))
}
