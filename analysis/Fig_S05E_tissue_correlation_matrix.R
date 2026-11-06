# Supplementary Fig. 5e | Tissue-level Pearson correlations
# Calculate Pearson's correlation for each of the ten distinct pairs of tissue-level
# measurements.
# P values are nominal across the ten comparisons. The observation unit is a tissue.
#
# Input: put Source_Data.xlsx beside this script or in the repository's data/ folder.
# Numerical packages: digest, dplyr, readxl.
# Plotting packages: dplyr, GGally, ggplot2, gtable, rlang, scales.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$pearson
#   analysis_result$tables$statistical_tests  # Effect-size estimates and reporting metadata
#   analysis_result$output_dir   # Folder containing Fig_S05E_correlation_matrix_statistics.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Supplementary Figure 5e: tissue marker counts and day-7 retention correlation matrix.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  library(dplyr)
  # Delta is the supplied difference between tracer-group medians at day 7 (%ID/g).
  plot_df <- readxl::read_excel(workbook, sheet = "Input S5d-e", range = "A1:F8", col_types = c("text", "numeric", "numeric", "numeric", "numeric", "numeric"), .name_repair = "minimal") %>%
    rename(
      Tissue = tissue,
      CD122 = No_of_CD122_positive_cell_by_tissue,
      CD25  = No_of_CD25_positive_cell_by_tissue,
      CD3   = No_of_CD3_positive_cell_by_tissue,
      Foxp3 = No_of_Foxp3_positive_cell_by_tissue,
      Delta = conc_difference_at_d7
    ) %>%
    mutate(Tissue = gsub("_", " ", Tissue))

  plot_vars <- c("CD122", "CD25", "CD3", "Foxp3", "Delta")

  # Labels describe the tissue-level input columns without changing their values.
  variable_descriptions <- c(
    CD122 = "Number of CD122 (Il2rb)-positive atlas cells per tissue",
    CD25 = "Number of CD25-positive atlas cells per tissue",
    CD3 = "Number of CD3-positive atlas cells per tissue",
    Foxp3 = "Number of Foxp3-positive atlas cells per tissue",
    Delta = "Day-7 median SLC-3010 uptake minus median hnTCB2 uptake (%ID/g)"
  )

  # Use one row per tissue and the complete cases available for each pair.
  if (anyNA(plot_df$Tissue) || anyDuplicated(plot_df$Tissue))
    stop("Expected one record per tissue in the correlation matrix.")
  pairs <- combn(plot_vars, 2, simplify = FALSE)
  comparisons <- vector("list", length(pairs))
  statistical_rows <- vector("list", length(pairs))
  fits <- vector("list", length(pairs))
  for (i in seq_along(pairs)) {
    pair <- pairs[[i]]
    data <- data.frame(x = plot_df[[pair[1]]], y = plot_df[[pair[2]]])
    complete <- is.finite(data$x) & is.finite(data$y)
    data <- data[complete, , drop = FALSE]
    if (nrow(data) < 3L || stats::sd(data$x) == 0 || stats::sd(data$y) == 0)
      stop("Cannot estimate a tissue correlation for ", paste(pair, collapse = " and "), ".")
    test <- stats::cor.test(data$x, data$y, method = "pearson", alternative = "two.sided", conf.level = 0.95)
    comparisons[[i]] <- data.frame(x = pair[1], y = pair[2], n = nrow(data),
      pearson_r = unname(test$estimate), pearson_p = test$p.value)
    # Report the correlation coefficient as the effect-size estimate.
    # These intervals and t statistics are returned by the same Pearson test.
    has_ci <- length(test$conf.int) == 2L
    statistical_rows[[i]] <- data.frame(
      Figure = "Supplementary Figure 5", Panel = "S5e",
      Comparison = paste(pair, collapse = " versus "),
      Variable_X = unname(variable_descriptions[pair[1]]),
      Variable_Y = unname(variable_descriptions[pair[2]]),
      n_pairs = nrow(data), Sample_unit = "Tissue",
      Tissues = paste(as.character(plot_df$Tissue[complete]), collapse = "; "),
      Measurement_design = "One paired tissue-level observation per tissue; the same tissue types are reused across the ten correlations",
      Test = "Pearson product-moment correlation", Alternative = "two-sided", Correction = "none",
      Statistic_name = names(test$statistic)[1L], Statistic_value = unname(test$statistic),
      df = unname(test$parameter), Effect_size_name = "Pearson r", Effect_size = unname(test$estimate),
      Effect_unit = "dimensionless", CI_level = if (has_ci) attr(test$conf.int, "conf.level") else NA_real_,
      CI_lower = if (has_ci) test$conf.int[1L] else NA_real_,
      CI_upper = if (has_ci) test$conf.int[2L] else NA_real_,
      CI_method = if (has_ci) "Fisher z transformation; unadjusted confidence interval returned by stats::cor.test" else
        "Not returned by stats::cor.test for this sample size",
      P_value = test$p.value, Adjusted_P_value = NA_real_,
      P_value_method = "Student t reference distribution with n - 2 degrees of freedom",
      R_workflow = 'stats::cor.test(x, y, method = "pearson", alternative = "two.sided", conf.level = 0.95)',
      Source_workbook = "Source_Data.xlsx", Source_sheet_range = "'Input S5d-e'!A1:F8",
      Source_script = "analysis/Fig_S05E_tissue_correlation_matrix.R",
      R_version = as.character(getRversion()), stats_version = as.character(getNamespaceVersion("stats")),
      Notes = paste("The correlation coefficient is the reported effect-size estimate.",
                    "Complete finite pairs are used for each correlation.",
                    "P values and confidence intervals are unadjusted across the ten comparisons.",
                    "The plotted straight lines are separate descriptive linear fits."),
      stringsAsFactors = FALSE
    )

    fit <- stats::lm(y ~ x, data = data)
    fits[[i]] <- data.frame(x = pair[1], y = pair[2], n = nrow(data),
      intercept = unname(stats::coef(fit)[1]), slope = unname(stats::coef(fit)[2]))
  }
  stat_tbl <- dplyr::bind_rows(comparisons)
  statistical_tests <- dplyr::bind_rows(statistical_rows)
  utils::write.csv(statistical_tests, file.path(output_dir, "Fig_S05E_Statistical_tests.csv"),
                   row.names = FALSE, na = "")
  trend_lines <- dplyr::bind_rows(fits)
  utils::write.csv(plot_df, file.path(output_dir, "Fig_S05E_tissue_measurements.csv"), row.names = FALSE)
  utils::write.csv(stat_tbl, file.path(output_dir, "Fig_S05E_correlation_matrix_statistics.csv"), row.names = FALSE)
  utils::write.csv(trend_lines, file.path(output_dir, "Fig_S05E_descriptive_linear_fits.csv"), row.names = FALSE)
  list(tables = list(pearson = stat_tbl, tissue_measurements = plot_df, linear_fits = trend_lines,
                     statistical_tests = statistical_tests))
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result, display = FALSE) {
  library(ggplot2)
  font_family <- "sans"
  plot_df <- result$tables$tissue_measurements
  plot_vars <- c("CD122", "CD25", "CD3", "Foxp3", "Delta")
  stat_tbl <- result$tables$pearson
  trends <- result$tables$linear_fits
  plot_labels <- c(
    "# CD122\npositive cells",
    "# CD25\npositive cells",
    "# CD3\npositive cells",
    "# Foxp3\npositive cells",
    "Δ concentration\n(%ID/g, day 7)"
  )


  p_stars <- function(p) {
    dplyr::case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ ""
    )
  }

  p_label <- function(p) {
    ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
  }

  correlation_theme <- function(base_size = if (display) 10 else 6.0) {
    theme_minimal(base_family = font_family, base_size = base_size) +
      theme(
        panel.grid.major = element_line(color = "#E6E9EF", linewidth = 0.25),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(color = "#111111", size = if (display) 9 else 5.0),
        axis.ticks = element_line(color = "#111111", linewidth = 0.25),
        panel.border = element_rect(color = "#BDBDBD", linewidth = 0.30, fill = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
  }

  lower_lm_panel <- function(data, mapping, ...) {
    x <- rlang::as_name(mapping$x)
    y <- rlang::as_name(mapping$y)
    fit <- trends[trends$x == x & trends$y == y, , drop = FALSE]
    if (nrow(fit) != 1L) stop("Missing precomputed linear fit for ", x, " and ", y)
    segment <- data.frame(x_value = range(data[[x]], na.rm = TRUE))
    segment$y_value <- fit$intercept + fit$slope * segment$x_value
    ggplot(data = data, mapping = mapping) +
      geom_point(size = 1.05, shape = 16, color = "#111111", alpha = 0.95) +
      geom_line(data = segment, aes(x = x_value, y = y_value), inherit.aes = FALSE,
                color = "#E41A1C", linewidth = 0.38) +
      scale_x_continuous(expand = expansion(mult = 0.08), breaks = scales::breaks_extended(n = 3)) +
      scale_y_continuous(expand = expansion(mult = 0.08), breaks = scales::breaks_extended(n = 3)) +
      correlation_theme()
  }

  diag_density_panel <- function(data, mapping, ...) {
    x <- GGally::eval_data_col(data, mapping$x)
    ggplot(data.frame(x = x), aes(x = x)) +
      geom_density(fill = "#D9D9D9", alpha = 0.55,
                   color = "#111111", linewidth = 0.35, adjust = 1.15) +
      scale_x_continuous(expand = expansion(mult = 0.08), breaks = scales::breaks_extended(n = 3)) +
      correlation_theme() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }

  upper_corr_panel <- function(data, mapping, ...) {
    x <- rlang::as_name(mapping$x)
    y <- rlang::as_name(mapping$y)
    row <- stat_tbl[(stat_tbl$x == x & stat_tbl$y == y) | (stat_tbl$x == y & stat_tbl$y == x), , drop = FALSE]
    if (nrow(row) != 1L) stop("Missing precomputed Pearson correlation for ", x, " and ", y)
    r <- row$pearson_r
    p <- row$pearson_p
    star <- p_stars(p)

    bg <- scales::col_numeric(
      palette = c("#2166AC", "#F7F7F7", "#B2182B"),
      domain = c(-1, 1)
    )(r)

    ggplot() +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
               fill = scales::alpha(bg, 0.11), color = NA) +
      annotate("text", x = 0.5, y = 0.60,
               label = sprintf("r = %.3f%s", r, star),
               family = font_family, fontface = "bold", size = if (display) 3.4 else 2.0,
               color = ifelse(p < 0.05, "#E41A1C", "#3A3A3A")) +
      annotate("text", x = 0.5, y = 0.34,
               label = p_label(p),
               family = font_family, size = if (display) 3.0 else 1.55, color = "#666666") +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      theme_void() +
      theme(panel.border = element_rect(color = "#BDBDBD", linewidth = 0.30, fill = NA))
  }

  # Plot
  p <- GGally::ggpairs(
    plot_df,
    columns = plot_vars,
    columnLabels = plot_labels,
    upper = list(continuous = upper_corr_panel),
    lower = list(continuous = lower_lm_panel),
    diag  = list(continuous = diag_density_panel),
    axisLabels = "none",
    progress = FALSE
  ) +
    theme(
      strip.background = element_rect(fill = "#F2F2F2", color = "#7A7A7A", linewidth = 0.35),
      strip.text = element_text(family = font_family, face = "bold", size = if (display) 9.8 else 5.4, color = "#111111", lineheight = 0.90),
      panel.spacing = grid::unit(0.7, "mm"),
      plot.margin = margin(2, 2, 8, 2, "pt")
    )

  add_caption_to_ggmatrix <- function(ggm, caption, caption_height = grid::unit(if (display) 32 else 8, "pt")) {
    gt <- GGally::ggmatrix_gtable(ggm)
    cap <- grid::textGrob(
      if (display) paste(strwrap(caption, width = 110), collapse = "\n") else caption,
      x = grid::unit(0, "npc"), y = grid::unit(0.5, "npc"),
      hjust = 0, vjust = 0.5,
      gp = grid::gpar(fontfamily = font_family, fontsize = if (display) 9 else 5.2, col = "#3A3A3A")
    )
    gt <- gtable::gtable_add_rows(gt, heights = caption_height, pos = nrow(gt))
    gt <- gtable::gtable_add_grob(gt, cap, t = nrow(gt), l = 1, r = ncol(gt), clip = "off")
    gt
  }

  matrix_caption <- sprintf(
    "Pearson correlation across %d tissues. Red line: linear fit. P values are unadjusted. * p < 0.05, ** p < 0.01, *** p < 0.001.",
    nrow(plot_df))

  list("Figure S5e – tissue marker and retention correlations" = add_caption_to_ggmatrix(p, matrix_caption))
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
  cat("Supplementary Fig. 5e | Tissue-level Pearson correlations\nSeven tissues; P values are nominal across the ten variable pairs. Full statistical records and confidence intervals are saved separately.\n")
  show("pearson", "All ten pairwise Pearson correlations", c("x", "y", "n", "pearson_r", "pearson_p"))
  cat("\nComplete result tables: ", paste(names(result$tables), collapse = ", "), "\n", sep = "")
  cat("Results directory: ", result$output_dir, "\n", sep = "")
  invisible(result)
}

# Source calculates the results, prints the analysis summary and displays the plots.
if (isTRUE(getOption("SLC3010.autorun", TRUE))) {
  analysis_result <- figure_result <- NULL
  analysis_result <- local({
    needed <- c("digest", "dplyr", "readxl")
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
    paths <- c(frames[basename(frames) == "Fig_S05E_tissue_correlation_matrix.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_S05E_tissue_correlation_matrix.R"], join(getwd(), "Fig_S05E_tissue_correlation_matrix.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_S05E_tissue_correlation_matrix.R"))], 1L)
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
      run <- join(base, "Fig_S05E_tissue_correlation_matrix", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
      needed <- c("GGally", "dplyr", "ggplot2", "gtable", "rlang", "scales")
      missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
      if (length(missing)) stop("Install these plotting packages: ", paste(missing, collapse = ", "))
      if (interactive() && identical(Sys.getenv("RSTUDIO"), "1") && !identical(names(grDevices::dev.cur()), "RStudioGD"))
        grDevices::dev.new(noRStudioGD = FALSE)
      measure <- !interactive() && grDevices::dev.cur() == 1L
      if (measure) grDevices::pdf(NULL)
      plots <- tryCatch({
        saved_plots <- suppressPackageStartupMessages(plot_results(analysis_result))
        display_plots <- suppressPackageStartupMessages(plot_results(analysis_result, display = TRUE))
        saved_plots
      }, finally = if (measure) grDevices::dev.off())
      for (i in seq_along(plots)) {
        size <- analysis_result$sizes[[names(plots)[i]]]
        if (is.null(size)) size <- c(9, 8)
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
