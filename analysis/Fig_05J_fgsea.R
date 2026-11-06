# Fig. 5j | Reported pathway enrichment and independent comparison
# Read the stored pathway statistics, reconstruct the gene ranking, and calculate a separate
# fgsea comparison with the supplied GO gene sets.
# The plot uses the stored enrichment statistics. Recomputed estimates and their comparisons
# are written under reproducibility/.
#
# Input: put Source_Data_RNA_inputs.xlsx beside this script or in the repository's data/
# folder.
# Numerical packages: digest, dplyr, readr, readxl, stringr, tibble / fgsea, BiocParallel (Bioconductor).
# Plotting packages: dplyr, ggplot2, ggrepel, stringr.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$stored_pathways
#   analysis_result$output_dir   # Folder containing Fig_05J_results.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Figure 5j: full saved GO:BP fgsea result visualization.
# Inputs: RNA_fgsea_saved, RNA_fgsea_rank_input and frozen RNA_GO_reference in Source Data.
# Plots the stored pathway statistics and the prespecified pathway labels.
# The input compares Triple hIL-2 C1D4 versus wild-type C1D4.
# The figure plots the stored statistics. An independent rank reconstruction is compared with these results.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  # Excel ranges used for named datasets in this analysis.
  table_ranges <- c(
    "RNA_GO_reference" = "'GO_BP_reference'!A1:F7933",
    "RNA_fgsea_saved" = "'MP_CD8_GSEA_results'!A1:K318",
    "RNA_fgsea_rank_input" = "'MP_CD8_GSEA_rank_input'!A1:C460"
  )
  RESULTS_DATASET <- "RNA_fgsea_saved"

  SEED <- 20260703

  PADJ_CUTOFF <- 0.05

  NES_LINE <- 1

  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(stringr))
  OUT_DIR <- output_dir

  require_columns <- function(x, cols, filename) {
      missing <- setdiff(cols, names(x))
      if (length(missing)) {
          stop("Missing columns in ", filename, ": ", paste(missing, collapse = ", "), "\nExisting columns: ", paste(names(x), 
              collapse = ", "), call. = FALSE)
      }
  }

  flatten_lists <- function(x) {
      x <- as.data.frame(x)
      for (nm in names(x)) if (is.list(x[[nm]])) {
          x[[nm]] <- vapply(x[[nm]], function(z) paste(z, collapse = ";"), character(1))
      }
      x
  }

  write_table <- function(x, filename) {
      readr::write_csv(flatten_lists(x), file.path(OUT_DIR, filename), na = "")
  }

  clean_gobp_name <- function(x) str_to_lower(str_replace_all(str_remove(x, "^GOBP_"), "_", " "))

  safe_neglog10 <- function(p) -log10(pmax(p, .Machine$double.xmin))

  fgsea_res <- readxl::read_excel(workbook, range = table_ranges[[RESULTS_DATASET]], .name_repair = "minimal")

  # Reconstruct ranks and compare reestimated enrichment statistics with the stored table.
  compare_enrichment <- function(workbook, table_ranges, fgsea_res, output_dir) {
    # Keep the caller's RNG state, including the state with no seed yet.
    had_rng_state <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    previous_rng_state <- if (had_rng_state) get(".Random.seed", envir = .GlobalEnv) else NULL
    previous_rng_kind <- RNGkind()
    on.exit({
      do.call(RNGkind, as.list(previous_rng_kind))
      if (had_rng_state) {
        assign(".Random.seed", previous_rng_state, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(list = ".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    read_gene_set_chunks <- function(dataset_id) {
      x <- as.data.frame(readxl::read_excel(workbook, range = table_ranges[[dataset_id]], col_types = "text", .name_repair = "minimal"), stringsAsFactors = FALSE)
      required <- c("Term", "Term_order", "Chunk", "Member_count", "Genes")
      if (!all(required %in% names(x))) stop("Incomplete membership table: ", dataset_id, call. = FALSE)
      for (n in c("Term_order", "Chunk", "Member_count")) {
        x[[n]] <- as.integer(x[[n]])
        if (anyNA(x[[n]])) stop("Invalid membership index in ", dataset_id, call. = FALSE)
      }
      if (anyDuplicated(x[c("Term_order", "Chunk")])) stop("Duplicated membership chunk.", call. = FALSE)
      x <- x[order(x$Term_order, x$Chunk), , drop = FALSE]
      members <- lapply(seq_len(nrow(x)), function(i) {
        g <- if (x$Member_count[i] == 0L) character() else strsplit(x$Genes[i], "|", fixed = TRUE)[[1]]
        if (length(g) != x$Member_count[i] || anyNA(g) || any(!nzchar(g)))
          stop("Invalid membership count in ", dataset_id, call. = FALSE)
        g
      })
      list(table = x, members = members)
    }

    read_go_reference <- function() {
      z <- read_gene_set_chunks("RNA_GO_reference")
      x <- z$table
      if (!"Source_rows" %in% names(x)) stop("Missing GO source-row indices.", call. = FALSE)
      row_ids <- lapply(seq_len(nrow(x)), function(i) {
        r <- as.integer(strsplit(x$Source_rows[i], "|", fixed = TRUE)[[1]])
        if (length(r) != x$Member_count[i] || anyNA(r)) stop("Invalid GO source-row indices.", call. = FALSE)
        r
      })
      d <- data.frame(gene_symbol = unlist(z$members, use.names = FALSE),
                      gs_name = rep(x$Term, x$Member_count),
                      source_row = unlist(row_ids, use.names = FALSE), stringsAsFactors = FALSE)
      d <- d[order(d$source_row), , drop = FALSE]
      if (!identical(d$source_row, seq_len(nrow(d)))) stop("GO source-row coverage is incomplete.", call. = FALSE)
      rownames(d) <- NULL
      d[, c("gene_symbol", "gs_name"), drop = FALSE]
    }

    write_enrichment_comparison <- function(calculated, stored, output_dir, prefix) {
      calculated <- as.data.frame(calculated, stringsAsFactors = FALSE)
      if (is.list(calculated$leadingEdge))
        calculated$leadingEdge <- vapply(calculated$leadingEdge, paste, character(1), collapse = ";")
      write.csv(calculated, file.path(output_dir, paste0(prefix, "_reestimated_results.csv")), row.names = FALSE, na = "")
      old <- as.data.frame(stored, stringsAsFactors = FALSE)
      if (!"pathway" %in% names(old) && "ID" %in% names(old)) names(old)[names(old) == "ID"] <- "pathway"
      if (!"pval" %in% names(old) && "pvalue" %in% names(old)) names(old)[names(old) == "pvalue"] <- "pval"
      if (!"padj" %in% names(old) && "p.adjust" %in% names(old)) names(old)[names(old) == "p.adjust"] <- "padj"
      if (!"ES" %in% names(old) && "enrichmentScore" %in% names(old)) names(old)[names(old) == "enrichmentScore"] <- "ES"
      if (anyDuplicated(calculated$pathway) || anyDuplicated(old$pathway)) stop("Duplicate pathway identifiers.")
      ids <- union(old$pathway, calculated$pathway)
      a <- match(ids, old$pathway)
      b <- match(ids, calculated$pathway)
      comparison <- data.frame(pathway = ids, in_stored = !is.na(a), in_reestimated = !is.na(b), stringsAsFactors = FALSE)
      coverage <- data.frame(stored_pathways = nrow(old), reestimated_pathways = nrow(calculated),
        shared_pathways = sum(!is.na(a) & !is.na(b)), stored_only = sum(!is.na(a) & is.na(b)),
        reestimated_only = sum(is.na(a) & !is.na(b)))
      write.csv(coverage, file.path(output_dir, paste0(prefix, "_pathway_coverage.csv")), row.names = FALSE)
      for (field in c("ES", "NES", "pval", "padj")) {
        comparison[[paste0(field, "_stored")]] <- old[[field]][a]
        comparison[[paste0(field, "_reestimated")]] <- calculated[[field]][b]
        comparison[[paste0(field, "_difference")]] <- calculated[[field]][b] - old[[field]][a]
      }
      write.csv(comparison, file.path(output_dir, paste0(prefix, "_stored_vs_reestimated.csv")), row.names = FALSE, na = "")
      summary <- data.frame(field = c("ES", "NES", "pval", "padj"), matched_rows = 0L,
                            exact_matches = 0L, max_absolute_difference = NA_real_)
      for (i in seq_len(nrow(summary))) {
        field <- summary$field[i]
        x <- comparison[[paste0(field, "_difference")]]
        finite <- is.finite(x)
        summary$matched_rows[i] <- sum(finite)
        summary$exact_matches[i] <- sum(x[finite] == 0)
        if (any(finite)) summary$max_absolute_difference[i] <- max(abs(x[finite]))
      }
      write.csv(summary, file.path(output_dir, paste0(prefix, "_comparison_summary.csv")), row.names = FALSE, na = "")
      invisible(summary)
    }

    write_method_record <- function(output_dir, prefix, notes) {
      writeLines(c(notes, "", capture.output(utils::sessionInfo())),
                 file.path(output_dir, paste0(prefix, "_method_and_runtime.txt")), useBytes = TRUE)
    }

    comparison_output_dir <- file.path(output_dir, "reproducibility")
    dir.create(comparison_output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(comparison_output_dir))
      stop("Cannot create the reproducibility output directory.", call. = FALSE)

    suppressPackageStartupMessages(library(dplyr))

    d <- readxl::read_excel(workbook, sheet = "MP_CD8_GSEA_rank_input", range = "A1:C460", .name_repair = "minimal")
    rank_table <- d %>% mutate(Gene_Symbol = as.character(Gene_Symbol),
                               rank_stat = as.numeric(`Signed_fold_change`), raw_p = as.numeric(`Raw_P_value`)) %>%
      filter(!is.na(Gene_Symbol), Gene_Symbol != "", is.finite(rank_stat))
    if (any(rank_table$raw_p < 0 | rank_table$raw_p > 1, na.rm = TRUE)) stop("Invalid nominal P value.")
    rank_table <- rank_table %>% arrange(desc(abs(rank_stat))) %>% group_by(Gene_Symbol) %>%
      slice(1) %>% ungroup() %>% arrange(desc(rank_stat))
    ranks <- setNames(rank_table$rank_stat, rank_table$Gene_Symbol)
    ranks <- sort(ranks + seq_along(ranks) * 1e-12, decreasing = TRUE)
    reference <- read_go_reference() %>% select(gs_name, gene_symbol) %>% distinct()
    pathways <- split(reference$gene_symbol, reference$gs_name)
    comparison_backend <- BiocParallel::SerialParam()
    set.seed(20260703)
    write.csv(data.frame(index = seq_along(.Random.seed), value = .Random.seed),
              file.path(comparison_output_dir, "Fig05J_initial_rng_state.csv"), row.names = FALSE)
    reestimated <- fgsea::fgseaMultilevel(pathways = pathways, stats = ranks, minSize = 10, maxSize = 500, eps = 0,
                                         BPPARAM = comparison_backend)
    write.csv(data.frame(Gene_Symbol = names(ranks), rank_stat = as.numeric(ranks)),
              file.path(comparison_output_dir, "Fig05J_ranked_genes.csv"), row.names = FALSE)
    stored <- fgsea_res
    comparison <- write_enrichment_comparison(reestimated, stored, comparison_output_dir, "Fig05J")
    write_method_record(comparison_output_dir, "Fig05J", c(
      "Independent enrichment estimates from the supplied ranked genes and GO memberships.",
      "These estimates do not reproduce the supplied Figure 5j statistics and do not replace the plotted results.",
      "The signed-linear fold changes for Triple hIL-2 C1D4 versus wild-type C1D4 are used directly as the ranking metric.",
      "Recorded call: fgseaMultilevel(minSize=10, maxSize=500, eps=0, BPPARAM=BiocParallel::SerialParam()); set.seed(20260703).",
      paste0("Enrichment backend: ", class(comparison_backend)[1L], "; workers=", BiocParallel::bpworkers(comparison_backend), "."),
      "The package versions and settings recorded below describe these independent estimates, not the analysis that produced the supplied results."))
    comparison
  }
  comparison <- compare_enrichment(workbook, table_ranges, fgsea_res, output_dir)
  result_columns <- c("pathway", "term", "NES", "pval", "padj", "ES", "log2err", "size", "leadingEdge", "direction", "minus_log10_padj")
  require_columns(fgsea_res, result_columns, RESULTS_DATASET)
  fgsea_res <- fgsea_res[, result_columns]

  fgsea_tbl <- tibble::as_tibble(fgsea_res)

  if ("selection_reason" %in% names(fgsea_tbl)) {
      stop("A selected-terms dataset cannot replace the complete stored result dataset.", call. = FALSE)
  }

  require_columns(fgsea_tbl, c("pathway", "NES", "padj"), "fgsea results")

  fgsea_tbl <- fgsea_tbl %>% mutate(pathway = as.character(pathway), NES = suppressWarnings(as.numeric(NES)), padj = suppressWarnings(as.numeric(padj)))

  if (anyDuplicated(fgsea_tbl$pathway)) stop("Duplicate pathway IDs in results. Supply one complete fgsea analysis, not concatenated analyses.")

  if (any(fgsea_tbl$padj < 0 | fgsea_tbl$padj > 1, na.rm = TRUE)) stop("Adjusted P values must lie in [0,1].")

  fgsea_tbl <- fgsea_tbl %>% mutate(term = clean_gobp_name(pathway), minus_log10_padj = safe_neglog10(padj), direction = case_when(is.finite(padj) & 
      padj < PADJ_CUTOFF & NES > 0 ~ "sig-Activated", is.finite(padj) & padj < PADJ_CUTOFF & NES < 0 ~ "sig-Repressed", 
      TRUE ~ "Not significant")) %>% arrange(desc(NES), padj)

  write_table(fgsea_tbl, "Fig_05J_results.csv")
  finite <- is.finite(fgsea_tbl$NES) & is.finite(fgsea_tbl$padj) & !is.na(fgsea_tbl$pathway) & nzchar(fgsea_tbl$pathway)
  direction <- factor(fgsea_tbl$direction[finite], levels = c("sig-Activated", "sig-Repressed", "Not significant"))
  counts <- as.data.frame(table(direction), responseName = "Pathways")
  list(tables = list(direction_counts = counts, enrichment_comparison = comparison,
                     stored_pathways = fgsea_tbl),
       threshold = c(adjusted_p = PADJ_CUTOFF, nes = NES_LINE))
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result) {
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(dplyr))
  fgsea_tbl <- result$tables$stored_pathways
  PADJ_CUTOFF <- unname(result$threshold["adjusted_p"])
  NES_LINE <- unname(result$threshold["nes"])
  LABEL_TOP_ACTIVATED <- 7L

  LABEL_TOP_REPRESSED <- 3L

  PRIORITY_TERMS <- c("cell division", "chromosome segregation", "mitotic nuclear division", "nuclear division", 
      "organelle fission", "cell cycle process", "mitotic cell cycle process", "mitotic sister chromatid segregation", 
      "nuclear chromosome segregation")

  REQUIRED_LABEL_TERMS <- c("cell division", "cell cycle process", "mitotic cell cycle process", "organelle fission", 
      "mitotic nuclear division", "mitotic sister chromatid segregation", "nuclear chromosome segregation", "chromosome segregation", 
      "cell cycle", "mitotic cell cycle", "cytokine production", "inorganic ion transmembrane transport", "cation transmembrane transport")

  REQUIRE_ALL_LABELS <- TRUE
  valid <- is.finite(fgsea_tbl$NES) & is.finite(fgsea_tbl$padj) & !is.na(fgsea_tbl$pathway) & nzchar(fgsea_tbl$pathway)

  if (any(!valid)) {
      warning(sum(!valid), " pathways have non-finite NES/adjusted P and are omitted from graphics only.")
  }

  plot_tbl <- fgsea_tbl[valid, , drop = FALSE]

  if (!nrow(plot_tbl)) stop("No fgsea results with finite NES and adjusted P.")

  priority_df <- plot_tbl %>% filter(direction == "sig-Activated", term %in% PRIORITY_TERMS)

  top_activated_df <- plot_tbl %>% filter(direction == "sig-Activated") %>% arrange(desc(minus_log10_padj), desc(NES)) %>% 
      slice_head(n = LABEL_TOP_ACTIVATED)

  top_repressed_df <- plot_tbl %>% filter(direction == "sig-Repressed") %>% arrange(desc(minus_log10_padj), NES) %>% 
      slice_head(n = LABEL_TOP_REPRESSED)

  negative_fallback_used <- nrow(top_repressed_df) < LABEL_TOP_REPRESSED

  if (negative_fallback_used) {
      top_repressed_df <- plot_tbl %>% filter(NES <= -NES_LINE) %>% arrange(desc(minus_log10_padj), NES) %>% slice_head(n = LABEL_TOP_REPRESSED)
  }

  required_names <- unique(stringr::str_to_lower(trimws(REQUIRED_LABEL_TERMS)))

  missing_required <- setdiff(required_names, plot_tbl$term)

  if (REQUIRE_ALL_LABELS && length(missing_required)) {
      stop("Figure 5j label(s) missing from the complete pathway table: ", paste(missing_required, collapse = "; "), 
          ". Supply the matching complete stored result dataset.",
          call. = FALSE)
  }

  required_label_df <- plot_tbl %>% filter(term %in% required_names)

  label_df <- bind_rows(priority_df, top_activated_df, top_repressed_df, required_label_df) %>% distinct(pathway, 
      .keep_all = TRUE) %>% mutate(is_ns = direction == "Not significant", label_term = stringr::str_wrap(paste0(term, 
      if_else(is_ns, " (NS)", "")), width = 30))

  stopifnot(!anyDuplicated(label_df$pathway))

  if (REQUIRE_ALL_LABELS && !all(required_names %in% label_df$term)) {
      stop("A required Figure 5j pathway is missing from the label set.")
  }

  if (any(is.na(label_df$label_term) | !nzchar(label_df$label_term))) {
      stop("A selected pathway has an empty label.")
  }
  # The panel uses stored pathway statistics. The independent fgsea results are exported in reproducibility/.
  direction_levels <- c("sig-Activated", "sig-Repressed", "Not significant")
  direction_labels <- c("Activated", "Repressed", "Not significant")
  plot_tbl$direction <- factor(plot_tbl$direction, levels = direction_levels)
  pathway_plot <- ggplot(plot_tbl, aes(minus_log10_padj, NES)) +
      geom_hline(yintercept = c(-NES_LINE, NES_LINE), linewidth = 0.3, linetype = "22", colour = "grey70") +
      geom_vline(xintercept = -log10(PADJ_CUTOFF), linewidth = 0.3, colour = "grey70") +
      geom_point(data = subset(plot_tbl, direction == "Not significant"), aes(colour = direction), size = 1, alpha = 0.65) +
      geom_point(data = subset(plot_tbl, direction != "Not significant"), aes(colour = direction), size = 2,
                 alpha = 0.90, show.legend = TRUE) +
      ggrepel::geom_text_repel(data = label_df, aes(label = label_term), size = 3, seed = 3010,
          min.segment.length = 0, colour = "grey15", max.overlaps = Inf, box.padding = 0.4, show.legend = FALSE) +
      scale_colour_manual("Pathway", values = c("sig-Activated"="#DF5B5F", "sig-Repressed"="#64B9D6",
          "Not significant"="#BEBEBE"), breaks = direction_levels, labels = direction_labels, drop = FALSE) +
      labs(x = expression(-log[10]~"adjusted P"), y = "Normalized enrichment score",
           caption = "Plotted: stored GO:BP pathway results. Independently recomputed estimates are reported separately.") +
      theme_classic(base_size = 10) + theme(legend.position = "bottom", plot.caption = element_text(hjust = 0))
  list("Figure 5j | Stored MP-CD8 pathway enrichment" = pathway_plot)
}

# Print the analysis estimates and tests without truncating to the first few observations.
print_analysis_results <- function(result) {
  old <- options(width = max(120L, getOption("width")))
  on.exit(options(old), add = TRUE)
  show <- function(key, label, columns) {
    x <- result$tables[[key]]
    if (!is.data.frame(x)) stop("Result table not available: ", key, call. = FALSE)
    missing <- setdiff(columns, names(x))
    if (length(missing)) stop("Missing report columns in ", key, ": ", paste(missing, collapse = ", "), call. = FALSE)
    x <- as.data.frame(x[, columns, drop = FALSE])
    schedule <- result$tables$treatment_schedule
    if ("group" %in% names(x) && is.data.frame(schedule) &&
        all(c("group", "first_admin", "second_admin") %in% names(schedule))) {
      x$group <- paste(schedule$first_admin, schedule$second_admin, sep = " -> ")[match(x$group, schedule$group)]
      names(x)[names(x) == "group"] <- "administration_sequence"
    }
    if (identical(key, "representative_pathways") && "pvalue" %in% names(x))
      x$pvalue <- ifelse(is.finite(x$pvalue) & x$pvalue <= 1e-10, "<= 1e-10", format(x$pvalue, digits = 6, scientific = TRUE))
    p_columns <- names(x)[grepl("(^p$|(^|[._])p([._]|$)|pvalue)", names(x), ignore.case = TRUE)]
    for (name in p_columns)
      if (is.numeric(x[[name]])) x[[name]] <- format.pval(x[[name]], digits = 6, eps = 0)
    cat("\n", label, "\n", sep = "")
    if (!nrow(x)) cat("No eligible rows.\n") else print(x, row.names = FALSE, digits = 6, max = .Machine$integer.max)
  }
  cat("Fig. 5j | Reported pathway enrichment and independent comparison\nThe figure uses the supplied pathway statistics. Independently recomputed statistics are separate reproducibility outputs and do not replace them.\n")
  show("direction_counts", "Pathway counts by reported enrichment direction", c("direction", "Pathways"))
  show("enrichment_comparison", "Independent calculation versus stored statistics", c("field", "matched_rows", "exact_matches", "max_absolute_difference"))
  cat("\nComplete result tables: ", paste(names(result$tables), collapse = ", "), "\n", sep = "")
  cat("Results directory: ", result$output_dir, "\n", sep = "")
  invisible(result)
}

# Source calculates the results, prints the analysis summary and displays the plots.
if (isTRUE(getOption("SLC3010.autorun", TRUE))) {
  analysis_result <- figure_result <- NULL
  analysis_result <- local({
    needed <- c("BiocParallel", "digest", "dplyr", "fgsea", "readr", "readxl", "stringr", "tibble")
    missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) stop("Install these analysis packages: ", paste(missing, collapse = ", "), call. = FALSE)
    join <- function(base, ...) file.path(sub("/+$", "", chartr("\\", "/", base)), ...)
    frames <- as.character(unlist(lapply(rev(sys.frames()), function(f) {
      name <- f[["ofile"]]
      if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) return(NULL)
      old <- f[["owd"]]
      if (length(old) == 1L && !is.na(old) && nzchar(old) && !grepl("^(/|[A-Za-z]:|\\\\|~)", name)) c(join(old, name), name) else name
    }), use.names = FALSE))
    cli <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
    direct <- !length(frames) && length(cli) && !interactive()
    editor <- if (interactive() && identical(Sys.getenv("RSTUDIO"), "1") && requireNamespace("rstudioapi", quietly = TRUE))
      tryCatch(as.character(rstudioapi::getSourceEditorContext()$path), error = function(e) character())
    paths <- c(frames[basename(frames) == "Fig_05J_fgsea.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_05J_fgsea.R"], join(getwd(), "Fig_05J_fgsea.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_05J_fgsea.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply one workbook path or directory, or use Source in RStudio.", call. = FALSE)
    path <- path.expand(if (length(args)) args[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR"))
    book <- "Source_Data_RNA_inputs.xlsx"
    sheets <- c("GO_BP_reference", "MP_CD8_GSEA_rank_input", "MP_CD8_GSEA_results")
    expected <- "2779beb9064ccd6090031bb707590c4de1d32a53e7b71326388a806ec7e9fd26"
    checksum <- function(f) digest::digest(readBin(f, "raw", n = file.info(f)$size), algo = "sha256", serialize = FALSE)
    if (!nzchar(path)) {
      candidates <- unique(join(c(here, join(repo, "data"), getwd(), join(getwd(), "data")), book))
      candidates <- candidates[file.exists(candidates) & !dir.exists(candidates)]
      rank <- vapply(candidates, function(f) tryCatch({
        if (all(sheets %in% readxl::excel_sheets(f))) 1L + identical(checksum(f), expected) else 0L
      }, error = function(e) 0L), integer(1))
      if (any(rank > 0L)) path <- candidates[which.max(rank)] else if (interactive()) {
        message("Select ", book); path <- file.choose()
      } else stop("Place ", book, " in data/ or beside this script, or pass its path.", call. = FALSE)
    }
    if (dir.exists(path)) path <- join(path, book)
    if (!file.exists(path) || dir.exists(path)) stop("Workbook not found: ", path, call. = FALSE)
    path <- normalizePath(path, winslash = "/", mustWork = TRUE)
    absent <- setdiff(sheets, readxl::excel_sheets(path))
    if (length(absent)) stop("Missing worksheets: ", paste(absent, collapse = ", "), call. = FALSE)
    sha <- checksum(path)
    status <- if (is.na(expected)) "Input workbook" else if (identical(sha, expected)) "Official journal workbook" else "Exploratory input (different workbook checksum)"
    is_root <- function(f) grepl("^/?$|^[A-Za-z]:/?$|^//[^/]+(/[^/]+)?/?$", chartr("\\", "/", f))
    preferred <- if (length(repo) && !is_root(dirname(repo))) join(dirname(repo), paste0(basename(repo), "_results")) else
      if (!length(repo) && !is_root(here)) join(here, "results")
    explicit <- Sys.getenv("SLC_OUTPUT_DIR")
    bases <- if (nzchar(explicit)) path.expand(explicit) else c(preferred, tools::R_user_dir("SLC3010", "data"))
    out <- NULL
    for (base in unique(bases)) {
      run <- join(base, "Fig_05J_fgsea", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
      if (!dir.create(run, recursive = TRUE, showWarnings = FALSE)) next
      probe <- file.path(run, ".write-test")
      if (file.create(probe, showWarnings = FALSE)) {
        unlink(probe)
        out <- normalizePath(run, winslash = "/")
        break
      }
    }
    if (is.null(out)) stop("Cannot write results under ", paste(bases, collapse = " or "), ". Set SLC_OUTPUT_DIR.", call. = FALSE)
    message(status, ": ", path, "\nResults directory: ", out)
    result <- suppressPackageStartupMessages(analyze_data(path, out))
    c(result, list(input = list(path = path, sha256 = sha, status = status), output_dir = out, tables_only = tables_only))
  })
  report_lines <- capture.output(print_analysis_results(analysis_result))
  cat(paste(report_lines, collapse = "\n"), "\n", sep = "")
  writeLines(report_lines, file.path(analysis_result$output_dir, "analysis_summary.txt"), useBytes = TRUE)
  if (!analysis_result$tables_only) {
    analysis_result$plots <- tryCatch(local({
      needed <- c("dplyr", "ggplot2", "ggrepel", "stringr")
      missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
      if (length(missing)) stop("Install these plotting packages: ", paste(missing, collapse = ", "))
      if (interactive() && identical(Sys.getenv("RSTUDIO"), "1") && !identical(names(grDevices::dev.cur()), "RStudioGD"))
        grDevices::dev.new(noRStudioGD = FALSE)
      measure <- !interactive() && grDevices::dev.cur() == 1L
      if (measure) grDevices::pdf(NULL)
      plots <- tryCatch(suppressPackageStartupMessages(plot_results(analysis_result)), finally = if (measure) grDevices::dev.off())
      for (i in seq_along(plots)) {
        size <- analysis_result$sizes[[names(plots)[i]]]
        if (is.null(size)) size <- c(9, 5.5)
        p <- plots[[i]]
        if (is.function(p)) p <- grid::grid.grabExpr(p(), width = size[1L], height = size[2L])
        if (interactive() && grid::is.grob(p)) { grid::grid.newpage(); grid::grid.draw(p) } else if (interactive()) print(p)
        ggplot2::ggsave(file.path(analysis_result$output_dir, sprintf("plot_%02d.png", i)), p,
                        width = size[1L], height = size[2L], units = "in", dpi = 160, bg = "white")
      }
      plots
    }), error = function(e) stop("Analysis tables saved in ", analysis_result$output_dir, ". Plotting failed: ", conditionMessage(e), call. = FALSE))
    figure_result <- analysis_result
  }
}
