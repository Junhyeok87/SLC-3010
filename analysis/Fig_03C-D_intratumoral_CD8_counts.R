# Fig. 3c-d | Intratumoral transferred and host CD8 counts
# For each readout, compare individual-mouse counts per gram of tumor using two-sided,
# unpaired equal-variance t tests and report differences, 95% confidence intervals and
# unadjusted P values.
# The reported difference is Vehicle minus the indicated drug, in millions of cells per gram.
#
# Input: put Source_Data.xlsx beside this script or in the repository's data/ folder.
# Numerical packages: digest, dplyr, readxl, stringr.
# Plotting packages: dplyr, ggplot2, patchwork, tibble.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$comparisons
#   analysis_result$output_dir   # Folder containing Fig_03C-D_vehicle_ttest_stats.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Fig. 3c,d: Tumor-weight-normalized transferred and host CD8 T-cell counts.
# Statistics: unpaired two-sided Student's t-tests; equal variances; no multiplicity adjustment.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(stringr))

  # Read the transferred and host per-mouse CD8 measurements.

  transferred <- readxl::read_excel(workbook, sheet = "Figure 3", range = "A781:C796", .name_repair = "minimal")
  host <- readxl::read_excel(workbook, sheet = "Figure 3", range = "A802:C817", .name_repair = "minimal")
  raw_data <- dplyr::bind_rows(
      dplyr::mutate(transferred, readout = "Thy1.1+CD8"),
      dplyr::mutate(host, readout = "Thy1.1-CD8")) %>%
      dplyr::transmute(mouse_id = mouse_id, treatment = treatment, readout = readout,
          value = cell_count_1e6_per_g_tumor)

  required_cols <- c("mouse_id", "treatment", "readout", "value")

  missing_cols <- setdiff(required_cols, colnames(raw_data))

  if (length(missing_cols) > 0) {
      stop("Input file is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }

  # Validate group and readout labels and retain the within-group ordering.

  treatment_levels <- c("Vehicle", "Beta IL-2", "SLC-3010")

  readout_levels <- c("Thy1.1+CD8", "Thy1.1-CD8")

  measurements <- raw_data %>%
    mutate(mouse_id = as.character(mouse_id), treatment = str_trim(as.character(treatment)),
           treatment = str_replace_all(treatment, "\\s*\\+\\s*", "+"),
           treatment = if_else(treatment == "BetaIL-2", "Beta IL-2", treatment),
           readout = str_trim(as.character(readout)), value = as.numeric(value)) %>%
    filter(!is.na(treatment), !is.na(readout), !is.na(value))
  if (any(!measurements$treatment %in% treatment_levels) ||
      any(!measurements$readout %in% readout_levels)) stop("Unknown treatment or readout.")
  if (anyNA(measurements$mouse_id) || anyDuplicated(measurements[c("mouse_id", "readout")]))
    stop("Each mouse must have one measurement per readout.")
  measurements <- measurements %>%
    mutate(treatment = factor(treatment, levels = treatment_levels),
           readout = factor(readout, levels = readout_levels)) %>%
    arrange(readout, treatment, suppressWarnings(as.numeric(mouse_id)), mouse_id)
  summary <- measurements %>%
    group_by(readout, treatment) %>%
    summarise(n = n(), mean = mean(value), median = median(value), sd = sd(value),
              sem = sd / sqrt(n), q1 = quantile(value, 0.25), q3 = quantile(value, 0.75),
              min = min(value), max = max(value), .groups = "drop")

  # The estimate, t statistic and 95% interval are defined in the same direction.
  compare_vehicle <- function(readout_id, active) {
    subset <- measurements[measurements$readout == readout_id, ]
    x <- subset$value[subset$treatment == "Vehicle"]
    y <- subset$value[subset$treatment == active]
    test <- tryCatch(stats::t.test(x, y, alternative = "two.sided", paired = FALSE,
                                  var.equal = TRUE, conf.level = 0.95), error = identity)
    ok <- !inherits(test, "error")
    data.frame(readout = readout_id, group1 = "Vehicle", group2 = active,
               n1 = length(x), n2 = length(y), mean1 = mean(x), mean2 = mean(y),
               effect_direction = "group1 minus group2", mean_difference = mean(x) - mean(y),
               confidence_level = 0.95, ci_low = if (ok) test$conf.int[1L] else NA_real_,
               ci_high = if (ok) test$conf.int[2L] else NA_real_,
               t_statistic = if (ok) unname(test$statistic) else NA_real_,
               df = if (ok) unname(test$parameter) else NA_real_,
               p_value = if (ok) test$p.value else NA_real_,
               test_method = "Student's t-test", alternative = "two.sided", paired = FALSE,
               var_equal = TRUE, p_adjustment = "none",
               test_status = if (ok) "Estimated" else conditionMessage(test))
  }
  comparisons <- do.call(rbind, lapply(readout_levels, function(r)
    do.call(rbind, lapply(treatment_levels[-1L], function(g) compare_vehicle(r, g)))))
  write.csv(measurements, file.path(output_dir, "Fig_03C-D_individual_measurements.csv"), row.names = FALSE)
  write.csv(summary, file.path(output_dir, "Fig_03C-D_summary_boxplot_stats_tumor_weight.csv"), row.names = FALSE)
  write.csv(comparisons, file.path(output_dir, "Fig_03C-D_vehicle_ttest_stats.csv"), row.names = FALSE)
  list(tables = list(comparisons = comparisons, summary = summary, observations = measurements))
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result) {
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(dplyr))

  observations <- result$tables$observations
  treatment_levels <- levels(observations$treatment)
  treatment_display_labels <- treatment_levels
  treatment_axis <- tibble::tibble(treatment = treatment_levels, treatment_label = treatment_levels,
                                   treatment_order = seq_along(treatment_levels), treatment_x = seq_along(treatment_levels))
  readout_levels <- levels(observations$readout)
  plot_data <- mutate(observations, treatment_x = as.numeric(treatment))

  point_jitter_width <- 0.085

  plot_data <- plot_data %>%  
      group_by(readout, treatment) %>% mutate(replicate_index = row_number(), n_in_group = n(), jitter_x = if_else(n_in_group <= 
      1, 0, (replicate_index - (n_in_group + 1)/2) * (2 * point_jitter_width/pmax(n_in_group - 1, 1))), x_point = treatment_x + 
      jitter_x) %>% ungroup()

  stat_comparisons <- tibble::tribble(~readout, ~group1, ~group2, ~y_top, ~y_drop, ~label_y, "Thy1.1+CD8", "Vehicle", 
      "Beta IL-2", 0.114, 0.012, 0.138, "Thy1.1+CD8", "Vehicle", "SLC-3010", 0.547, 0.012, 0.571, "Thy1.1-CD8", "Vehicle", 
      "Beta IL-2", 17.8, 0.85, 19.5, "Thy1.1-CD8", "Vehicle", "SLC-3010", 32.4, 0.85, 34.1)

  p_to_asterisks <- function(p) dplyr::case_when(is.na(p) ~ "NA", p < 1e-04 ~ "****", p < 0.001 ~ "***", p < 0.005 ~ 
      "**", p < 0.05 ~ "*", TRUE ~ "ns")

  stat_data <- left_join(stat_comparisons, result$tables$comparisons,
                         by = c("readout", "group1", "group2")) %>%
    mutate(label = p_to_asterisks(p_value),
           x1 = match(group1, treatment_levels), x2 = match(group2, treatment_levels),
           x_mid = (x1 + x2) / 2, y_low = y_top - y_drop)

  base_family <- "Arial"
  base_size_single <- 6.8
  show_panel_title <- TRUE
  panel_title_size <- 7.2

  readout_display_labels <- c(`Thy1.1+CD8` = "Thy1.1⁺CD8⁺", `Thy1.1-CD8` = "Thy1.1⁻CD8⁺")

  y_axis_titles <- c(`Thy1.1+CD8` = "Thy1.1⁺CD8⁺ (×10⁶)\n/ Tumor weight (g)", `Thy1.1-CD8` = "Thy1.1⁻CD8⁺ (×10⁶)\n/ Tumor weight (g)")

  y_axis_limits <- list(`Thy1.1+CD8` = c(0, 0.6), `Thy1.1-CD8` = c(0, 40))

  y_axis_breaks <- list(`Thy1.1+CD8` = seq(0, 0.6, by = 0.2), `Thy1.1-CD8` = seq(0, 40, by = 10))

  treatment_colors <- c(Vehicle = "#4D4D4D", `Beta IL-2` = "#3B5FD9", `SLC-3010` = "#E53935")

  treatment_fills <- c(Vehicle = "white", `Beta IL-2` = "#3B5FD9", `SLC-3010` = "#E53935")

  treatment_shapes <- c(Vehicle = 21, `Beta IL-2` = 17, `SLC-3010` = 15)

  point_size_single <- 1.45
  point_stroke <- 0.44
  box_width <- 0.5
  box_alpha <- 0.16

  theme_fig3cd_single <- function(base_size = base_size_single, base_family = base_family) theme_classic(base_size = base_size, 
      base_family = base_family) %+replace% theme(strip.background = element_blank(), strip.text = element_blank(), 
      plot.title.position = "panel", plot.title = element_text(size = panel_title_size, face = "plain", color = "black", 
          hjust = 0.5, lineheight = 1, margin = margin(b = 1.1, unit = "pt")), axis.title.x = element_blank(), axis.title.y = element_text(size = base_size + 
          0.2, face = "plain", color = "black", angle = 90, margin = margin(r = 2.2)), axis.text.x = element_text(size = base_size, 
          color = "black", angle = 45, hjust = 1, vjust = 1, margin = margin(t = 2.1)), axis.text.y = element_text(size = base_size, 
          color = "black", margin = margin(r = 1.2)), axis.line = element_line(color = "black", linewidth = 0.34), 
      axis.ticks = element_line(color = "black", linewidth = 0.34), axis.ticks.length = unit(1.4, "pt"), panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), legend.position = "none", plot.margin = margin(t = 0.7, r = 0.8, b = 5, 
          l = 0.6, unit = "pt"))

  add_common_scales <- function(p) p + scale_x_continuous(breaks = treatment_axis$treatment_x, labels = treatment_axis$treatment_label, 
      limits = c(0.5, length(treatment_levels) + 0.5), expand = expansion(mult = c(0, 0))) + scale_color_manual(values = treatment_colors, 
      breaks = treatment_levels, labels = treatment_display_labels, drop = FALSE) + scale_fill_manual(values = treatment_fills, 
      breaks = treatment_levels, labels = treatment_display_labels, guide = "none", drop = FALSE) + scale_shape_manual(values = treatment_shapes, 
      breaks = treatment_levels, labels = treatment_display_labels, guide = "none", drop = FALSE) + coord_cartesian(clip = "off")

  make_single_panel_plot <- function(readout_id, point_size = point_size_single) {
      d <- plot_data %>% filter(as.character(readout) == readout_id)
      p <- ggplot(d) + geom_boxplot(aes(x = treatment_x, y = value, group = treatment, fill = treatment), orientation = "x", 
          width = box_width, outlier.shape = NA, color = "black", linewidth = 0.3, alpha = box_alpha) + geom_point(aes(x = x_point, 
          y = value, color = treatment, fill = treatment, shape = treatment), size = point_size, stroke = point_stroke, 
          alpha = 1) + geom_segment(data = stat_data %>% filter(readout == readout_id), aes(x = x1, xend = x2, y = y_top, 
          yend = y_top), inherit.aes = FALSE, linewidth = 0.28, color = "black") + geom_segment(data = stat_data %>% 
          filter(readout == readout_id), aes(x = x1, xend = x1, y = y_low, yend = y_top), inherit.aes = FALSE, linewidth = 0.28, 
          color = "black") + geom_segment(data = stat_data %>% filter(readout == readout_id), aes(x = x2, xend = x2, 
          y = y_low, yend = y_top), inherit.aes = FALSE, linewidth = 0.28, color = "black") + geom_text(data = stat_data %>% 
          filter(readout == readout_id), aes(x = x_mid, y = label_y, label = label), inherit.aes = FALSE, family = base_family, 
          fontface = "bold", size = 2.2, lineheight = 0.95, color = "black") + scale_y_continuous(limits = y_axis_limits[[readout_id]], 
          breaks = y_axis_breaks[[readout_id]], expand = expansion(mult = c(0, 0.045))) + labs(x = NULL, y = y_axis_titles[[readout_id]], 
          title = if (isTRUE(show_panel_title)) 
              unname(readout_display_labels[[readout_id]])
          else NULL)
      add_common_scales(p) + theme_fig3cd_single(base_size = base_size_single, base_family = base_family)
  }

  fig3c_vertical_boxplot <- make_single_panel_plot("Thy1.1+CD8", point_size = point_size_single)

  fig3d_vertical_boxplot <- make_single_panel_plot("Thy1.1-CD8", point_size = point_size_single)

  panels <- patchwork::wrap_plots(fig3c_vertical_boxplot, fig3d_vertical_boxplot, nrow = 1)

  list("Fig. 3c–d | Intratumoral CD8 counts" = panels)
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
  cat("Fig. 3c-d | Intratumoral transferred and host CD8 counts\nUnpaired, two-sided Student t-tests with equal variances; P values are unadjusted. Difference: Vehicle minus drug.\n")
  show("comparisons", "Mean differences and 95% confidence intervals (millions/g)", c("readout", "group1", "group2", "n1", "n2", "mean_difference", "ci_low", "ci_high"))
  show("comparisons", "Student t-tests", c("readout", "group1", "group2", "t_statistic", "df", "p_value"))
  show("summary", "Group means and variability", c("readout", "treatment", "n", "mean", "sd", "sem"))
  cat("\nComplete result tables: ", paste(names(result$tables), collapse = ", "), "\n", sep = "")
  cat("Results directory: ", result$output_dir, "\n", sep = "")
  invisible(result)
}

# Source calculates the results, prints the analysis summary and displays the plots.
if (isTRUE(getOption("SLC3010.autorun", TRUE))) {
  analysis_result <- figure_result <- NULL
  analysis_result <- local({
    needed <- c("digest", "dplyr", "readxl", "stringr")
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
    paths <- c(frames[basename(frames) == "Fig_03C-D_intratumoral_CD8_counts.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_03C-D_intratumoral_CD8_counts.R"], join(getwd(), "Fig_03C-D_intratumoral_CD8_counts.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_03C-D_intratumoral_CD8_counts.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply one workbook path or directory, or use Source in RStudio.", call. = FALSE)
    path <- path.expand(if (length(args)) args[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR"))
    book <- "Source_Data.xlsx"
    sheets <- c("Figure 3")
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
      run <- join(base, "Fig_03C-D_intratumoral_CD8_counts", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
      needed <- c("dplyr", "ggplot2", "patchwork", "tibble")
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
