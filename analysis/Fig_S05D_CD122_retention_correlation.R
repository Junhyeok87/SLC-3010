# Supplementary Fig. 5d | Tissue CD122 and tracer-retention correlation
# Across seven tissues, calculate Spearman's correlation between Il2rb-positive cell counts
# and median SLC-3010 uptake minus median hnTCB2 uptake.
# Each point represents a tissue. The straight line is a descriptive fit, separate from the
# rank correlation.
#
# Input: put Source_Data.xlsx beside this script or in the repository's data/ folder.
# Numerical packages: digest, dplyr, readr, readxl, stringr, tibble.
# Plotting packages: dplyr, ggplot2, ggrepel, scales.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$spearman
#   analysis_result$tables$statistical_tests  # Effect-size estimate and reporting metadata
#   analysis_result$output_dir   # Folder containing Fig_S05D_CD122_correlation_statistics.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Supplementary Figure 5d: tissue CD122-positive cell counts and day-7 retention.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  library(dplyr)
  # The study input supplies day-7 median(SLC-3010) minus median(hnTCB2), in %ID/g.
  pretty_tissue <- function(x) {
    x %>%
      as.character() %>%
      stringr::str_replace_all("_", " ") %>%
      stringr::str_replace_all("Large Intestine", "Large intestine") %>%
      stringr::str_replace_all("large intestine", "Large intestine") %>%
      stringr::str_replace_all("Urinary bladder", "Bladder") %>%
      stringr::str_replace_all("UB", "Bladder")
  }
  raw <- readxl::read_excel(workbook, sheet = "Input S5d-e", range = "A1:F8", col_types = c("text", "numeric", "numeric", "numeric", "numeric", "numeric"), .name_repair = "minimal")

  corr <- raw %>%
    transmute(
      tissue = pretty_tissue(.data$tissue),
      cd122_pos_cells = as.numeric(.data$No_of_CD122_positive_cell_by_tissue),
      conc_difference_at_d7 = as.numeric(.data$conc_difference_at_d7)
    ) %>%
    filter(!is.na(tissue), !is.na(cd122_pos_cells), !is.na(conc_difference_at_d7)) %>%
    mutate(
      tissue = factor(
        tissue,
        levels = c("Spleen", "Lung", "Heart", "Liver", "Large intestine", "Bladder", "Kidney")
      )
    )

  if (anyNA(corr$tissue) || anyDuplicated(corr$tissue))
    stop("Expected one row for each recognized tissue in the correlation table.")
  # Spearman statistics are computed from ranks.
  sp <- suppressWarnings(stats::cor.test(corr$cd122_pos_cells, corr$conc_difference_at_d7,
                                         method = "spearman", alternative = "two.sided", exact = FALSE))

  stats_tbl <- tibble::tibble(
    method = "Spearman",
    estimate = unname(sp$estimate),
    p_value = sp$p.value,
    n = nrow(corr)
  )

  # Report the rank correlation coefficient as the effect-size estimate.
  # n is the number of paired tissue-level observations, not a count of mice or cells.
  statistical_tests <- data.frame(
    Figure = "Supplementary Figure 5", Panel = "S5d",
    Comparison = "CD122-positive atlas cell count versus day-7 tracer-retention difference",
    Variable_X = "Number of CD122 (Il2rb)-positive atlas cells per tissue",
    Variable_Y = "Day-7 median SLC-3010 uptake minus median hnTCB2 uptake (%ID/g)",
    n_pairs = nrow(corr), Sample_unit = "Tissue",
    Tissues = paste(as.character(corr$tissue), collapse = "; "),
    Measurement_design = "One paired tissue-level observation per tissue; the tissue types are also used in S5e",
    Test = "Spearman rank correlation", Alternative = "two-sided", Correction = "none",
    Statistic_name = names(sp$statistic)[1L], Statistic_value = unname(sp$statistic),
    df = NA_real_, Effect_size_name = "Spearman rho", Effect_size = unname(sp$estimate),
    Effect_unit = "dimensionless", CI_level = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_,
    CI_method = "Not returned by stats::cor.test for the Spearman method; no additional interval was calculated",
    P_value = sp$p.value, Adjusted_P_value = NA_real_,
    P_value_method = "Approximate Spearman P value (exact = FALSE)",
    R_workflow = 'stats::cor.test(x, y, method = "spearman", alternative = "two.sided", exact = FALSE)',
    Source_workbook = "Source_Data.xlsx", Source_sheet_range = "'Input S5d-e'!A1:F8",
    Source_script = "analysis/Fig_S05D_CD122_retention_correlation.R",
    R_version = as.character(getRversion()), stats_version = as.character(getNamespaceVersion("stats")),
    Notes = paste("The correlation coefficient is the reported effect-size estimate.",
                  "The S statistic and P value are returned by stats::cor.test.",
                  "Degrees of freedom are not returned for this method.",
                  "The plotted straight line is a separate descriptive linear fit; no multiplicity adjustment is applied."),
    stringsAsFactors = FALSE
  )
  utils::write.csv(statistical_tests, file.path(output_dir, "Fig_S05D_Statistical_tests.csv"),
                   row.names = FALSE, na = "")

  readr::write_csv(corr, file.path(output_dir, "Fig_S05D_tissue_retention.csv"))
  readr::write_csv(stats_tbl, file.path(output_dir, "Fig_S05D_CD122_correlation_statistics.csv"))
  # The descriptive straight line is separate from the Spearman test.
  model <- stats::lm(conc_difference_at_d7 ~ cd122_pos_cells, data = corr)
  trend <- data.frame(n = nrow(corr), intercept = unname(stats::coef(model)[1]),
                      slope = unname(stats::coef(model)[2]))
  readr::write_csv(trend, file.path(output_dir, "Fig_S05D_descriptive_linear_fit.csv"))
  list(tables = list(spearman = stats_tbl, tissue_measurements = corr, linear_fit = trend,
                     statistical_tests = statistical_tests))
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result, display = FALSE) {
  library(ggplot2)
  library(dplyr)
  base_family <- "sans"
  corr <- result$tables$tissue_measurements
  stats_tbl <- result$tables$spearman
  linear_fit <- result$tables$linear_fit
  endpoints <- data.frame(x = range(corr$cd122_pos_cells))
  endpoints$y <- linear_fit$intercept + linear_fit$slope * endpoints$x
  p_format <- function(p) {
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  }
  base_theme <- function(base_size = 6.1, base_family = "sans") {
    theme_classic(base_size = base_size, base_family = base_family) +
      theme(
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major = element_line(colour = "#E7E8EF", linewidth = 0.23),
        panel.grid.minor = element_blank(),
        axis.line  = element_line(colour = "#222222", linewidth = 0.34),
        axis.ticks = element_line(colour = "#222222", linewidth = 0.30),
        axis.ticks.length = unit(1.3, "pt"),
        axis.text = element_text(colour = "#222222", size = 5.5),
        axis.title.x = element_text(size = 6.2, margin = margin(t = 4.0)),
        axis.title.y = element_text(size = 6.2, margin = margin(r = 4.0), lineheight = 0.88),
        legend.position = "none",
        plot.margin = margin(2.5, 5.0, 2.5, 4.0, "pt")
      )
  }

  point_size <- 1.45
  line_width <- 0.55
  label_size <- if (display) 3.4 else 2.25
  stat_size  <- if (display) 3.2 else 1.75
  stat_label <- sprintf("Spearman r = %.3f\np = %s", stats_tbl$estimate, p_format(stats_tbl$p_value))

  x_min <- min(corr$cd122_pos_cells, na.rm = TRUE)
  x_max <- max(corr$cd122_pos_cells, na.rm = TRUE)
  y_min <- min(corr$conc_difference_at_d7, na.rm = TRUE)
  y_max <- max(corr$conc_difference_at_d7, na.rm = TRUE)
  x_span <- x_max - x_min
  y_span <- y_max - y_min

  # Leave room for labels and the statistics.
  x_limits <- c(x_min - 0.05 * x_span, x_max + 0.10 * x_span)
  y_limits <- c(max(0, y_min - 0.14 * y_span), y_max + 0.30 * y_span)

  stat_x_main <- x_limits[2] - 0.035 * diff(x_limits)
  stat_y_main <- y_limits[1] + 0.115 * diff(y_limits)

  label_dat <- corr %>%
    mutate(
      nudge_x = dplyr::case_when(
        tissue == "Kidney" ~ 13,
        tissue == "Bladder" ~ 11,
        tissue == "Large intestine" ~ 18,
        tissue == "Heart" ~ 9,
        tissue == "Liver" ~ 7,
        tissue == "Lung" ~ 12,
        tissue == "Spleen" ~ -16,
        TRUE ~ 0
      ),
      nudge_y = dplyr::case_when(
        tissue == "Kidney" ~ 0.38,
        tissue == "Bladder" ~ -0.25,
        tissue == "Large intestine" ~ -0.36,
        tissue == "Heart" ~ 0.38,
        tissue == "Liver" ~ 0.32,
        tissue == "Lung" ~ -0.43,
        tissue == "Spleen" ~ -0.32,
        TRUE ~ 0
      )
    ) %>%
    arrange(tissue)

  p_main <- ggplot(corr, aes(x = cd122_pos_cells, y = conc_difference_at_d7)) +
    geom_line(data = endpoints, aes(x = x, y = y), inherit.aes = FALSE,
              linewidth = line_width, colour = "#2F62FF") +
    geom_point(size = point_size, colour = "#111111") +
    ggrepel::geom_text_repel(
      data = label_dat,
      aes(label = tissue),
      nudge_x = label_dat$nudge_x,
      nudge_y = label_dat$nudge_y,
      size = label_size,
      family = base_family,
      colour = "#111111",
      seed = 3010,
      min.segment.length = 0,
      segment.size = 0.20,
      segment.colour = "#444444",
      box.padding = 0.20,
      point.padding = 0.18,
      force = 3.0,
      force_pull = 0.55,
      max.iter = 50000,
      max.overlaps = Inf
    ) +
    annotate(
      "text",
      x = stat_x_main,
      y = stat_y_main,
      label = stat_label,
      hjust = 1,
      vjust = 0,
      size = stat_size,
      family = base_family,
      colour = "#111111"
    ) +
    scale_x_continuous(
      limits = x_limits,
      breaks = scales::breaks_pretty(n = 4),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = scales::breaks_pretty(n = 4),
      expand = expansion(mult = 0)
    ) +
    labs(
      x = expression("No. of CD122"^"+"~"cells"),
      y = "Concentration difference\nat day 7 (%ID/g)",
      caption = "Each point is one tissue. Retention: difference in day-7 tracer medians. Blue line: linear fit; statistic: Spearman correlation."
    ) +
    coord_cartesian(clip = "off") +
    base_theme(base_family = base_family)
  list("Figure S5d – tissue CD122 and day-7 retention" = p_main)
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
  cat("Supplementary Fig. 5d | Tissue CD122 and tracer-retention correlation\nThe unit of observation is a tissue. The straight line in the figure is a separate descriptive fit.\n")
  show("spearman", "Spearman correlation across tissues", c("method", "estimate", "p_value", "n"))
  show("tissue_measurements", "Seven tissue-level observations", c("tissue", "cd122_pos_cells", "conc_difference_at_d7"))
  cat("\nComplete result tables: ", paste(names(result$tables), collapse = ", "), "\n", sep = "")
  cat("Results directory: ", result$output_dir, "\n", sep = "")
  invisible(result)
}

# Source calculates the results, prints the analysis summary and displays the plots.
if (isTRUE(getOption("SLC3010.autorun", TRUE))) {
  analysis_result <- figure_result <- NULL
  analysis_result <- local({
    needed <- c("digest", "dplyr", "readr", "readxl", "stringr", "tibble")
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
    paths <- c(frames[basename(frames) == "Fig_S05D_CD122_retention_correlation.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_S05D_CD122_retention_correlation.R"], join(getwd(), "Fig_S05D_CD122_retention_correlation.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_S05D_CD122_retention_correlation.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply one workbook path or directory, or use Source in RStudio.", call. = FALSE)
    path <- path.expand(if (length(args)) args[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR"))
    book <- "Source_Data.xlsx"
    sheets <- c("Input S5d-e")
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
      run <- join(base, "Fig_S05D_CD122_retention_correlation", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
      needed <- c("dplyr", "ggplot2", "ggrepel", "scales")
      missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
      if (length(missing)) stop("Install these plotting packages: ", paste(missing, collapse = ", "))
      if (interactive() && identical(Sys.getenv("RSTUDIO"), "1") && !identical(names(grDevices::dev.cur()), "RStudioGD"))
        grDevices::dev.new(noRStudioGD = FALSE)
      measure <- !interactive() && grDevices::dev.cur() == 1L
      if (measure) grDevices::pdf(NULL)
      plots <- tryCatch({
        saved_plots <- suppressPackageStartupMessages(plot_results(analysis_result))
        display_plots <- lapply(suppressPackageStartupMessages(plot_results(analysis_result, display = TRUE)), make_display_plot)
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
