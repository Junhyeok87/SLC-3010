# Supplementary Fig. 5c | Il2rb-positive atlas cells by tissue
# Count Il2rb-positive cells by tissue and verify that each cell's stored expression equals
# log(count + 1).
# The worksheet contains positive cells. Their counts do not estimate the fraction of all
# atlas cells expressing Il2rb.
#
# Input: put Source_Data.xlsx beside this script or in the repository's data/ folder.
# Numerical packages: digest, dplyr, readr, readxl.
# Plotting packages: forcats, ggplot2, scales.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$tissue_counts
#   analysis_result$output_dir   # Folder containing Fig_S05C_CD122_positive_cells_by_tissue_summary.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Supplementary Figure 5c: nonzero Il2rb expression in atlas cells by tissue.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  library(dplyr)
  final_data <- readxl::read_excel(workbook, sheet = "Supplementary Figure 5", range = "A782:G1202", col_types = c("text", "text", "text", "text", "numeric", "numeric", "numeric"), .name_repair = "minimal")

  required_columns <- c("tissue", "cell_ontology_class", "T_cell", "mouse.id", "counts", "ln_(count+1)", "n")
  if (!all(required_columns %in% names(final_data))) stop("The saved atlas table has missing columns.")
  if (any(!is.finite(final_data$counts) | final_data$counts <= 0)) stop("Expected positive finite Il2rb counts.")
  if (!isTRUE(all.equal(final_data$`ln_(count+1)`, log(final_data$counts + 1), tolerance = 1e-12))) {
    stop("Stored transformed values do not match log(count + 1).")
  }
  final_data <- final_data %>%
    mutate(tissue = factor(tissue),
           cell_ontology_class = factor(cell_ontology_class),
           mouse.id = factor(mouse.id),
           T_cell = factor(T_cell, levels = c("T cell", "others")))
  count_final <- final_data %>% count(tissue, name = "n")
  count_final <- dplyr::arrange(count_final, dplyr::desc(n))
  cell_expression <- final_data
  cell_expression$n <- count_final$n[match(cell_expression$tissue, count_final$tissue)]
  readr::write_csv(cell_expression, file.path(output_dir, "Fig_S05C_Il2rb_positive_cells.csv"))
  readr::write_csv(count_final, file.path(output_dir, "Fig_S05C_CD122_positive_cells_by_tissue_summary.csv"))
  list(tables = list(tissue_counts = count_final, cell_expression = cell_expression))
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result) {
  library(ggplot2)
  final_data <- result$tables$cell_expression
  final_data$tissue_plot <- forcats::fct_reorder(final_data$tissue, final_data$n, .desc = TRUE)
  BASE_SIZE      <- 6.6
  AXIS_TEXT_SIZE <- 5.6
  AXIS_TITLE_SIZE <- 6.8
  LEGEND_TEXT_SIZE <- 6.5
  POINT_SIZE     <- 1.00
  POINT_ALPHA    <- 0.60
  JITTER_WIDTH   <- 0.20
  # Figure
  p1 <- ggplot(final_data, aes(x = tissue_plot, y = `ln_(count+1)`, color = T_cell)) +
    # Fix display-only jitter across image formats and repeated runs.
    geom_point(
      position = position_jitter(width = JITTER_WIDTH, height = 0, seed = 3010),
      alpha = POINT_ALPHA,
      size = POINT_SIZE,
      stroke = 0
    ) +
    scale_color_manual(
      values = c("T cell" = "red", "others" = "gray30"),
      labels = c("T cells", "Other cells"),
      name = ""
    ) +
    scale_x_discrete(labels = function(x) gsub("_", " ", x, fixed = TRUE)) +
    expand_limits(y = c(0, 10)) +
    scale_y_continuous(
      breaks = scales::breaks_width(2),
      expand = expansion(mult = c(0.02, 0.03))
    ) +
    labs(
      x = "Tissue",
      y = "CD122 expression, ln(count + 1)",
      caption = "Each point is one atlas cell with detectable Il2rb expression."
    ) +
    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(size = 1.8, alpha = 1)
      )
    ) +
    theme_bw(base_size = BASE_SIZE, base_family = "sans") +
    theme(
      legend.position = "top",
      legend.justification = "right",
      legend.box.just = "right",
      legend.direction = "horizontal",
      legend.box.margin = margin(0, 0, -2, 0),
      legend.margin = margin(0, 0, 0, 0),
      legend.key.width = grid::unit(0.16, "in"),
      legend.key.height = grid::unit(0.10, "in"),
      legend.text = element_text(size = LEGEND_TEXT_SIZE, color = "black"),

      axis.title.x = element_text(size = AXIS_TITLE_SIZE, color = "black", margin = margin(t = 1.5)),
      axis.title.y = element_text(size = AXIS_TITLE_SIZE, color = "black", margin = margin(r = 2.5)),
      axis.text.x = element_text(size = AXIS_TEXT_SIZE, color = "gray25", angle = 0, hjust = 0.5, vjust = 0.5),
      axis.text.y = element_text(size = AXIS_TEXT_SIZE, color = "gray25"),
      axis.ticks = element_line(linewidth = 0.25, color = "black"),
      axis.ticks.length = unit(1.3, "pt"),

      panel.border = element_rect(color = "black", linewidth = 0.35, fill = NA),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1, 2, 1, 1, "pt")
    )
  list("Figure S5c – CD122 expression by tissue" = p1)
}

# Larger screen text; the original plot remains the publication export.
make_display_plot <- function(p) {
  if (!inherits(p, "ggplot")) return(p)
  if (inherits(p, "patchwork"))
    p$patches$plots <- lapply(p$patches$plots, make_display_plot)
  screen_theme <- ggplot2::theme(
    text = ggplot2::element_text(size = 11),
    axis.title = ggplot2::element_text(size = 12),
    axis.title.x = ggplot2::element_text(size = 12, margin = ggplot2::margin(t = 6)),
    axis.title.y = ggplot2::element_text(size = 12, margin = ggplot2::margin(r = 6)),
    axis.text = ggplot2::element_text(size = 10),
    axis.text.x = ggplot2::element_text(size = 10),
    axis.text.y = ggplot2::element_text(size = 10),
    legend.title = ggplot2::element_text(size = 11),
    legend.text = ggplot2::element_text(size = 10),
    strip.text = ggplot2::element_text(size = 11),
    plot.title = ggplot2::element_text(size = 13),
    plot.subtitle = ggplot2::element_text(size = 11),
    plot.caption = ggplot2::element_text(size = 9),
    plot.margin = ggplot2::margin(7, 7, 7, 7))
  # Keep intentionally hidden axes, titles and guides hidden.
  for (name in names(screen_theme))
    if (inherits(p$theme[[name]], "element_blank")) screen_theme[[name]] <- p$theme[[name]]
  p <- p + screen_theme

  p
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
  cat("Supplementary Fig. 5c | Il2rb-positive atlas cells by tissue\nCounts refer to positive cells in the supplied Tabula Muris FACS subset, not the fraction of all atlas cells.\n")
  show("tissue_counts", "Numbers of Il2rb-positive cells", c("tissue", "n"))
  cat("\nComplete result tables: ", paste(names(result$tables), collapse = ", "), "\n", sep = "")
  cat("Results directory: ", result$output_dir, "\n", sep = "")
  invisible(result)
}

# Source calculates the results, prints the analysis summary and displays the plots.
if (isTRUE(getOption("SLC3010.autorun", TRUE))) {
  analysis_result <- figure_result <- NULL
  analysis_result <- local({
    needed <- c("digest", "dplyr", "readr", "readxl")
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
    paths <- c(frames[basename(frames) == "Fig_S05C_CD122_expression.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_S05C_CD122_expression.R"], join(getwd(), "Fig_S05C_CD122_expression.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_S05C_CD122_expression.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply one workbook path or directory, or use Source in RStudio.", call. = FALSE)
    path <- path.expand(if (length(args)) args[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR"))
    book <- "Source_Data.xlsx"
    sheets <- c("Supplementary Figure 5")
    expected <- "1f0e88ae7ed898b04b89cbcaa13249c0e70b76612b20df6c77e12b396cb2aa9c"
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
      run <- join(base, "Fig_S05C_CD122_expression", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
      needed <- c("forcats", "ggplot2", "scales")
      missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
      if (length(missing)) stop("Install these plotting packages: ", paste(missing, collapse = ", "))
      if (interactive() && identical(Sys.getenv("RSTUDIO"), "1") && !identical(names(grDevices::dev.cur()), "RStudioGD"))
        grDevices::dev.new(noRStudioGD = FALSE)
      measure <- !interactive() && grDevices::dev.cur() == 1L
      if (measure) grDevices::pdf(NULL)
      plots <- tryCatch({
        saved_plots <- suppressPackageStartupMessages(plot_results(analysis_result))
        display_plots <- lapply(saved_plots, make_display_plot)
        saved_plots
      }, finally = if (measure) grDevices::dev.off())
      for (i in seq_along(plots)) {
        size <- analysis_result$sizes[[names(plots)[i]]]
        if (is.null(size)) size <- c(9, 5.5)
        p <- plots[[i]]
        if (is.function(p)) p <- grid::grid.grabExpr(p(), width = size[1L], height = size[2L])
        display_plot <- display_plots[[i]]
        if (is.function(display_plot)) display_plot <- p
        if (interactive() && grid::is.grob(display_plot)) { grid::grid.newpage(); grid::grid.draw(display_plot) } else if (interactive()) print(display_plot)
        ggplot2::ggsave(file.path(analysis_result$output_dir, sprintf("plot_%02d.png", i)), p,
                        width = size[1L], height = size[2L], units = "in", dpi = 160, bg = "white")
      }
      structure(plots, display_plots = display_plots)
    }), error = function(e) stop("Analysis tables saved in ", analysis_result$output_dir, ". Plotting failed: ", conditionMessage(e), call. = FALSE))
    analysis_result$display_plots <- attr(analysis_result$plots, "display_plots")
    attr(analysis_result$plots, "display_plots") <- NULL
    figure_result <- analysis_result
  }
}
