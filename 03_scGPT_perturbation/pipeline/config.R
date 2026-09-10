# Shared R configuration loader; paths are relative to the YAML file.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || args[[1]] == "--help") {
  cat("Usage: Rscript SCRIPT.R config.yaml [--query]
")
  quit(status = if (length(args)) 0L else 1L)
}
config_path <- normalizePath(args[[1]], mustWork = TRUE)
cfg <- yaml::read_yaml(config_path)
resolve_path <- function(x) {
  if (grepl("^/|^[A-Za-z]:", x)) x else file.path(dirname(config_path), x)
}
