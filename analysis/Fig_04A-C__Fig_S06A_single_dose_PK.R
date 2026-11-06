# Fig. 4a-c and Supplementary Fig. 6a | Single-dose serum pharmacokinetics
# Calculate mean concentrations by dose, analyte and time, then obtain Cmax and trapezoidal
# AUC from each group mean profile over its observed interval.
# Positive observations define the first and last integration times, which can differ by
# analyte and dose.
#
# Input: put Source_Data.xlsx beside this script or in the repository's data/ folder.
# Numerical packages: digest, dplyr, readr, readxl.
# Plotting packages: dplyr, ggplot2, scales.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$exposure_metrics
#   analysis_result$output_dir   # Folder containing pk_metrics_from_mean_profiles.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Fig. 4a-c and Supplementary Fig. 6a: Single-dose serum profiles, mean-profile AUC and Cmax.
# Only positive, nonmissing observed concentrations enter the source profile means.
# AUC uses sorted linear trapezoids through available positive mean observations.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  PK_raw <- readxl::read_excel(workbook, sheet = "Figure 4", range = "A7:E300", .name_repair = "minimal")

  required_cols <- c("ID", "Time", "Dose", "Detection", "Conc")

  missing_cols <- setdiff(required_cols, names(PK_raw))

  if (length(missing_cols) > 0) {
      stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  PK <- dplyr::mutate(PK_raw, ID = factor(ID), Time = as.numeric(Time), Dose_num = as.numeric(Dose), Dose = factor(Dose_num, 
      levels = c(12, 36, 72), ordered = TRUE), Detection = factor(Detection, levels = c("SLC-3010", "hnTCB2")), Conc = as.numeric(Conc))

  PK_obs <- dplyr::filter(PK, !is.na(ID), !is.na(Time), !is.na(Dose), !is.na(Detection), !is.na(Conc), Conc > 0)

  summary_pk <- dplyr::summarise(dplyr::group_by(PK_obs, Detection, Dose, Dose_num, Time), n_obs = dplyr::n(), mean_conc = mean(Conc, 
      na.rm = TRUE), sd_conc = if (dplyr::n() > 1) stats::sd(Conc, na.rm = TRUE) else NA_real_, .groups = "drop")

  trap_auc <- function(time, conc) {
      keep <- is.finite(time) & is.finite(conc)
      time <- time[keep]
      conc <- conc[keep]
      ord <- order(time)
      time <- time[ord]
      conc <- conc[ord]
      if (length(time) < 2) 
          return(NA_real_)
      sum(diff(time) * (head(conc, -1) + tail(conc, -1))/2)
  }

  # Area under the observed group-mean profile in ng*h/mL.
  pk_metrics <- summary_pk |>
      dplyr::arrange(Detection, Dose_num, Time) |>
      dplyr::group_by(Detection, Dose, Dose_num) |>
      dplyr::summarise(Cmax = max(mean_conc), AUC_observed = trap_auc(Time, mean_conc),
          first_observed_h = min(Time), last_observed_h = max(Time), observed_timepoints = dplyr::n(),
          first_time_n = dplyr::first(n_obs), last_time_n = dplyr::last(n_obs), .groups = "drop")

  readr::write_csv(summary_pk, file.path(output_dir, "summary_concentration_by_group.csv"))

  readr::write_csv(pk_metrics, file.path(output_dir, "pk_metrics_from_mean_profiles.csv"))

  write.csv(summary_pk, file.path(output_dir, "Fig_04A-B_concentration_summary.csv"), row.names = FALSE)

  write.csv(pk_metrics, file.path(output_dir, "Fig_04C_Fig_S06A_PK_metrics.csv"), row.names = FALSE)

  list(tables = list(exposure_metrics = pk_metrics, concentration_summary = summary_pk),
       measurements = PK_obs)
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result) {
  base_family <- "Arial"
  PK_obs <- result$measurements
  summary_pk <- result$tables$concentration_summary
  pk_metrics <- result$tables$exposure_metrics

  x_breaks <- c(0, 24, 48, 72, 168, 336)

  x_lim <- c(0, 360)

  y_breaks <- c(10, 100, 1000, 10000, 1e+05)

  y_lim <- c(6, 180000)

  pk_cols <- c(`12` = "#111111", `36` = "#2867B2", `72` = "#C22D2D")

  metric_cols <- c(`SLC-3010` = "#4D4D4D", hnTCB2 = "#C51B36")

  metric_shapes <- c(`SLC-3010` = 16, hnTCB2 = 17)

  paper_theme <- function(base_size = 7.2) ggplot2::theme_classic(base_size = base_size, base_family = base_family) + 
      ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.line = ggplot2::element_line(linewidth = 0.38, colour = "black"), 
          axis.ticks = ggplot2::element_line(linewidth = 0.38, colour = "black"), axis.ticks.length = grid::unit(1.35, 
              "mm"), axis.text = ggplot2::element_text(size = base_size, colour = "black"), axis.title = ggplot2::element_text(size = base_size + 
              0.4, colour = "black"), plot.title = ggplot2::element_text(size = base_size + 0.8, face = "bold", hjust = 0.5, 
              margin = ggplot2::margin(b = 3)), legend.title = ggplot2::element_text(size = base_size), legend.text = ggplot2::element_text(size = base_size - 
              0.2), legend.key.height = grid::unit(3, "mm"), legend.key.width = grid::unit(4.2, "mm"), plot.margin = ggplot2::margin(4, 
              5, 5, 6, "pt"))

  make_pk_panel <- function(det, show_dose_legend = FALSE) {
      dat_obs <- dplyr::filter(PK_obs, Detection == det)
      dat_mean <- dplyr::filter(summary_pk, Detection == det)
      p <- ggplot2::ggplot() + ggplot2::geom_line(data = dat_obs, mapping = ggplot2::aes(x = Time, y = Conc, group = interaction(Dose, 
          ID), colour = Dose), linewidth = 0.22, alpha = 0.14, show.legend = FALSE) + ggplot2::geom_line(data = dat_mean, 
          mapping = ggplot2::aes(x = Time, y = mean_conc, group = Dose, colour = Dose), linewidth = 0.68, lineend = "round") + 
          ggplot2::geom_point(data = dat_mean, mapping = ggplot2::aes(x = Time, y = mean_conc, colour = Dose), size = 1.45) + 
          ggplot2::scale_x_continuous(limits = x_lim, breaks = x_breaks, expand = ggplot2::expansion(mult = c(0, 
              0.02))) + ggplot2::scale_y_log10(limits = y_lim, breaks = y_breaks, labels = scales::label_comma(accuracy = 1)) + 
          ggplot2::scale_colour_manual(name = "Dose (µg)", values = pk_cols, breaks = c("12", "36", "72"), labels = c("12", 
              "36", "72")) + ggplot2::labs(title = det, x = "Time (hours) post SLC-3010 injection", y = "Concentration (ng/mL)") + 
          ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1, byrow = TRUE, override.aes = list(linewidth = 0.75, 
              size = 1.6))) + paper_theme(7.2) + ggplot2::theme(axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)), 
          legend.position = if (show_dose_legend) 
              c(0.98, 0.98)
          else "none", legend.justification = c(1, 1), legend.direction = "horizontal", legend.background = ggplot2::element_rect(fill = scales::alpha("white", 
              0.88), colour = NA), legend.box.margin = ggplot2::margin(0, 0, 0, 0))
      p
  }

  make_metric_panel <- function(metric, title_expr, y_expr, breaks, divisor = 1, show_legend = FALSE, legend_inside = FALSE) {
      p <- ggplot2::ggplot(pk_metrics, ggplot2::aes(x = Dose_num, y = .data[[metric]], group = Detection, colour = Detection, 
          shape = Detection)) + ggplot2::geom_line(linewidth = 0.55, lineend = "round") + ggplot2::geom_point(size = 1.75) + 
          ggplot2::scale_x_continuous(breaks = c(12, 36, 72), limits = c(7, 77), expand = ggplot2::expansion(mult = c(0, 
              0))) + ggplot2::scale_y_continuous(breaks = breaks, labels = function(x) scales::number(x/divisor, 
          accuracy = 1), expand = ggplot2::expansion(mult = c(0.05, 0.08))) + ggplot2::scale_colour_manual(values = metric_cols) + 
          ggplot2::scale_shape_manual(values = metric_shapes) + ggplot2::labs(title = title_expr, x = "SLC-3010 dose (µg)", 
          y = y_expr, colour = NULL, shape = NULL) + ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1, byrow = TRUE, 
          override.aes = list(linewidth = 0.75, size = 1.8)), shape = ggplot2::guide_legend(nrow = 1, byrow = TRUE)) + 
          paper_theme(7.2) + ggplot2::theme(axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)), 
          legend.position = if (!show_legend) 
              "none"
          else if (legend_inside) 
              c(0.98, 0.98)
          else "top", legend.justification = if (legend_inside) 
              c(1, 1)
          else "center", legend.direction = "horizontal", legend.background = ggplot2::element_rect(fill = scales::alpha("white", 
              0.88), colour = NA), legend.box.margin = ggplot2::margin(0, 0, 0, 0))
      p
  }

  p_slc <- make_pk_panel("SLC-3010", show_dose_legend = TRUE)

  p_tcb <- make_pk_panel("hnTCB2", show_dose_legend = TRUE) + ggplot2::labs(title = "hnTCB2")

  p_auc <- make_metric_panel("AUC_observed", expression(AUC["observed"]), expression(paste(AUC["observed"], " (", 10^6, 
      " ng·h/mL)")), c(0, 1e+06, 2e+06, 3e+06, 4e+06), divisor = 1e+06, show_legend = TRUE, legend_inside = FALSE)

  p_cmax <- make_metric_panel("Cmax", NULL, expression(paste(C[max], " (", 10^3, " ng/mL)")), c(20000, 40000, 60000, 
      80000), divisor = 1000, show_legend = TRUE, legend_inside = TRUE) + ggplot2::theme(legend.text = ggplot2::element_text(size = 6.2), 
      legend.key.width = grid::unit(3.6, "mm"), legend.key.height = grid::unit(2.6, "mm"), legend.position = c(0.02, 
          0.98), legend.justification = c(0, 1), legend.direction = "vertical") + ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2, 
      byrow = TRUE, override.aes = list(linewidth = 0.75, size = 1.8)), shape = ggplot2::guide_legend(nrow = 2, byrow = TRUE))

  list("Fig. 4a - serum SLC-3010" = p_slc,
                    "Fig. 4b - serum hnTCB2" = p_tcb,
                    "Fig. 4c - mean-profile AUC" = p_auc,
                    "Supplementary Fig. 6a - mean-profile Cmax" = p_cmax)
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
  cat("Fig. 4a-c and Supplementary Fig. 6a | Single-dose serum pharmacokinetics\nCmax is in ng/mL and AUC in ng-hour/mL. AUC uses each mean profile's observed interval, which may differ between analytes and doses.\n")
  show("exposure_metrics", "Peak concentration and observed exposure from group mean profiles", c("Detection", "Dose", "Cmax", "AUC_observed", "first_observed_h", "last_observed_h", "observed_timepoints", "first_time_n", "last_time_n"))
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
    paths <- c(frames[basename(frames) == "Fig_04A-C__Fig_S06A_single_dose_PK.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_04A-C__Fig_S06A_single_dose_PK.R"], join(getwd(), "Fig_04A-C__Fig_S06A_single_dose_PK.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_04A-C__Fig_S06A_single_dose_PK.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply one workbook path or directory, or use Source in RStudio.", call. = FALSE)
    path <- path.expand(if (length(args)) args[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR"))
    book <- "Source_Data.xlsx"
    sheets <- c("Figure 4")
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
      run <- join(base, "Fig_04A-C__Fig_S06A_single_dose_PK", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
      needed <- c("dplyr", "ggplot2", "scales")
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
