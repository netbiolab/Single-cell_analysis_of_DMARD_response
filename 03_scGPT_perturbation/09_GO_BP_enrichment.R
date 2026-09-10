#!/usr/bin/env Rscript
# Purpose: GO Biological Process enrichment and q-value-sorted barplots of SPC top 20.
# Input: SRC selected-gene CSV and archived EnrichR workbook, or explicit live --query.
# Output: GO_BP EnrichR XLSX, input summary CSV and top-10 term PDF/PNG.
# Preserves GO_Biological_Process_2021, overlap fraction and original q-value sorting.
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(forcats)
  library(enrichR)
})

script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(script_file)), "pipeline", "config.R"))
input_csv <- file.path(resolve_path(cfg$paths$output_dir), "SRC", "spc_euclidean_top20_selected_genes.csv")
output_dir <- file.path(resolve_path(cfg$paths$output_dir), "GO_BP")
selected_perturbations <- unlist(cfg$figure5$targets)
db_configs <- list(list(db = cfg$gobp$database, short = "GO_BP"))
top_n <- cfg$gobp$top_terms
q_cut <- cfg$gobp$q_reference_line
query_live <- "--query" %in% args
cache_file <- resolve_path(cfg$paths$enrichr_cache)
pick_column <- function(df, candidates) {
  hits <- intersect(candidates, names(df))
  if (!length(hits)) NA_character_ else hits[1]
}

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  stats::reorder(paste(x, within, sep = sep), by, FUN = fun, ...)
}

strip_reorder_suffix <- function(x, sep = "___") {
  sub(paste0(sep, ".*$"), "", x)
}

extract_top_genes <- function(row_df) {
  gene_cols <- grep("^top[0-9]+_gene$", names(row_df), value = TRUE)
  gene_cols <- gene_cols[order(as.integer(str_match(gene_cols, "^top([0-9]+)_gene$")[, 2]))]
  genes <- unlist(row_df[1, gene_cols], use.names = FALSE)
  genes <- as.character(genes)
  genes <- genes[!is.na(genes) & genes != ""]
  unique(genes)
}

build_plot_df <- function(enr_df, perturb_gene, db_short) {
  term_col <- pick_column(enr_df, c("Term", "term", "Description", "Pathway", "Term.Name"))
  q_col <- pick_column(enr_df, c("Adjusted.P.value", "Adjusted.P.Value", "q.value", "qvalue", "Adj.P", "FDR"))
  p_col <- pick_column(enr_df, c("P.value", "P.Value", "pvalue", "P"))
  ol_col <- pick_column(enr_df, c("Overlap", "overlap"))

  if (is.na(term_col)) {
    stop(sprintf("Term column not found for %s enrichment result.", perturb_gene))
  }
  if (is.na(q_col)) {
    if (is.na(p_col)) {
      stop(sprintf("Adjusted p-value or p-value column not found for %s enrichment result.", perturb_gene))
    }
    q_col <- p_col
  }

  enr_df %>%
    transmute(
      term = .data[[term_col]],
      qval = as.numeric(.data[[q_col]]),
      overlap_raw = if (!is.na(ol_col)) .data[[ol_col]] else NA
    ) %>%
    mutate(term = str_remove(term, "\\s*\\(.*\\)$")) %>%
    mutate(
      term = if (db_short == "Reactome") {
        str_remove(term, "\\s*R-HSA-[0-9]+$")
      } else {
        term
      }
    ) %>%
    mutate(
      overlap_pct = if (!all(is.na(overlap_raw))) {
        xy <- str_split_fixed(as.character(overlap_raw), "/", 2)
        suppressWarnings(as.numeric(xy[, 1]) / as.numeric(xy[, 2]))
      } else {
        NA_real_
      }
    ) %>%
    mutate(neglog10q = -log10(pmax(qval, .Machine$double.xmin))) %>%
    arrange(desc(neglog10q), desc(overlap_pct)) %>%
    slice_head(n = top_n) %>%
    mutate(
      perturb_gene = perturb_gene,
      term_within = reorder_within(term, neglog10q, perturb_gene)
    )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)
if (query_live) setEnrichrSite("Enrichr")

spc_df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
target_df <- spc_df %>% filter(perturb_gene %in% selected_perturbations)

if (nrow(target_df) != length(selected_perturbations)) {
  missing_genes <- setdiff(selected_perturbations, target_df$perturb_gene)
  stop(sprintf("Missing perturbation rows in input CSV: %s", paste(missing_genes, collapse = ", ")))
}

input_summary_df <- bind_rows(lapply(selected_perturbations, function(gene) {
  row_df <- target_df %>% filter(perturb_gene == gene)
  top_genes <- extract_top_genes(row_df)
  data.frame(
    perturb_gene = gene,
    input_genes = paste(top_genes, collapse = ","),
    n_input_genes = length(top_genes),
    stringsAsFactors = FALSE
  )
}))

for (cfg in db_configs) {
  db_name <- cfg$db
  db_short <- cfg$short

  workbook_path <- file.path(output_dir, paste0("spc_top20_selected_genes_", db_short, "_EnrichR_qvalue_sorted.xlsx"))
  summary_csv_path <- file.path(output_dir, paste0("spc_top20_selected_genes_", db_short, "_input_summary_qvalue_sorted.csv"))
  plot_pdf_path <- file.path(output_dir, paste0("spc_top20_selected_genes_", db_short, "_enrichr_barplot_top10_qvalue_sorted.pdf"))
  plot_png_path <- file.path(output_dir, paste0("spc_top20_selected_genes_", db_short, "_enrichr_barplot_top10_qvalue_sorted.png"))

  wb <- createWorkbook()
  plot_dfs <- list()

  for (gene in selected_perturbations) {
    row_df <- target_df %>% filter(perturb_gene == gene)
    top_genes <- extract_top_genes(row_df)

    if (query_live) {
      enrich_list <- enrichr(top_genes, db_name)
      enrich_df <- enrich_list[[1]]
    } else {
      if (!file.exists(cache_file)) stop("Provide paths.enrichr_cache, or use --query for a new live analysis")
      enrich_df <- read.xlsx(cache_file, sheet = substr(paste0(db_short, "_", gene), 1, 31))
    }

    sheet_name <- substr(paste0(db_short, "_", gene), 1, 31)
    addWorksheet(wb, sheet_name)
    writeDataTable(wb, sheet = sheet_name, x = enrich_df, rowNames = TRUE, colNames = TRUE)

    plot_dfs[[gene]] <- build_plot_df(enrich_df, gene, db_short)
  }

  saveWorkbook(wb, workbook_path, overwrite = TRUE)
  write.csv(input_summary_df, summary_csv_path, row.names = FALSE)

  plot_df_all <- bind_rows(plot_dfs) %>%
    mutate(perturb_gene = factor(perturb_gene, levels = selected_perturbations))

  p <- ggplot(plot_df_all, aes(x = neglog10q, y = term_within, fill = overlap_pct)) +
    geom_col(width = 0.8, color = NA) +
    geom_vline(xintercept = -log10(q_cut), linetype = "dashed", color = "red", linewidth = 0.6) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_y_discrete(labels = strip_reorder_suffix) +
    scale_fill_gradientn(
      colours = c("#6B8EC1", "#8FA0CC", "#B59DC0", "#CF8FAA", "#DD7E96", "#E36F85", "#D85E73"),
      na.value = "grey85",
      labels = percent_format(accuracy = 1),
      name = "overlap\npercentage"
    ) +
    labs(
      title = paste0(db_short, " Enrichment of SPC Top-20 Genes"),
      subtitle = "Top 10 enriched terms per perturbation",
      x = expression(-log[10](qvalue)),
      y = db_short
    ) +
    facet_wrap(~ perturb_gene, ncol = 1, scales = "free_y") +
    theme_light(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title.y = element_text(margin = margin(r = 8)),
      axis.title.x = element_text(margin = margin(t = 6)),
      legend.title = element_text(lineheight = 0.9),
      strip.text = element_text(face = "bold")
    )

  ggsave(plot_pdf_path, p, width = 9, height = 12)
  ggsave(plot_png_path, p, width = 9, height = 12, dpi = 300)

  message("Saved workbook: ", workbook_path)
  message("Saved summary CSV: ", summary_csv_path)
  message("Saved plot PDF: ", plot_pdf_path)
  message("Saved plot PNG: ", plot_png_path)
}
