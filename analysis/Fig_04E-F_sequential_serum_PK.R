# Fig. 4e,f | Serum concentrations after sequential administration
# Source Data: Figure 4, A373:G668; concentrations in micrograms/mL.
# Run this file with Source in RStudio. Place Source_Data.xlsx beside the script
# or in the repository's data folder; otherwise select it when prompted.
# Results include concentration means/SEM and descriptive per-mouse PK summaries.
# AUC uses the linear trapezoidal rule over 0, 0.5, 4, 8 and 24 hours.
# PK group summaries use complete profiles only, without imputing missing values.
# Cmax/Tmax describe the observed sampling grid; no model or group test is fitted.
# Numerical packages: digest, dplyr, readxl. Plotting packages: ggplot2, patchwork.

analyze_data <- function(workbook, output_dir) {
  suppressPackageStartupMessages(library(dplyr))
  group_levels <- c("G1", "G2", "G3", "G4", "G5")
  readout_levels <- c("SLC-3010", "hnTCB2")
  expected_times <- c(0, 0.5, 4, 8, 24)
  raw_data <- readxl::read_excel(workbook, sheet = "Figure 4", range = "A373:G668", .name_repair = "minimal")
  required <- c("mouse_id", "group", "first_admin", "second_admin", "readout", "time_h", "value")
  absent <- setdiff(required, names(raw_data))
  if (length(absent)) stop("Input is missing columns: ", paste(absent, collapse = ", "), call. = FALSE)
  raw_data$plot_label <- paste0(raw_data$first_admin, "\u2192", raw_data$second_admin)
  plot_data <- raw_data %>%
    select(mouse_id, group, first_admin, second_admin, plot_label, readout, time_h, value) %>%
    mutate(mouse_id = as.character(mouse_id), group = factor(as.character(group), levels = group_levels),
           first_admin = as.character(first_admin), second_admin = as.character(second_admin),
           plot_label = as.character(plot_label), readout = factor(as.character(readout), levels = readout_levels),
           time_h = as.numeric(time_h), value = as.numeric(value)) %>%
    filter(!is.na(group), !is.na(readout), !is.na(time_h), !is.na(value))
  if (anyDuplicated(plot_data[c("readout", "group", "mouse_id", "time_h")]))
    stop("Repeated readout/group/mouse/time records in the input.", call. = FALSE)
  if (any(!is.finite(plot_data$value)) || any(!plot_data$time_h %in% expected_times))
    stop("Expected finite concentrations at 0, 0.5, 4, 8 or 24 hours.", call. = FALSE)

  # Retain all observed measurements for the figure's mean and SEM.
  summary_data <- plot_data %>% group_by(readout, group, time_h) %>%
    summarise(n = n(), mean = mean(value, na.rm = TRUE),
              sd = ifelse(n > 1, sd(value, na.rm = TRUE), 0), sem = sd/sqrt(n),
              ymin = pmax(mean - sem, 0), ymax = mean + sem, .groups = "drop") %>%
    mutate(readout = factor(as.character(readout), levels = readout_levels),
           group = factor(as.character(group), levels = group_levels))
  treatments <- plot_data %>% distinct(group, first_admin, second_admin, plot_label)
  if (anyDuplicated(treatments$group)) stop("A group has inconsistent treatment labels.", call. = FALSE)
  concentration_summary <- summary_data %>% left_join(treatments, by = "group") %>%
    select(readout, group, first_admin, second_admin, treatment = plot_label, time_h, n, mean, sd, sem)

  # Incomplete profiles remain visible in the subject table but do not supply
  # AUC, Cmax or Tmax to the complete-profile group summaries.
  subject_pk <- plot_data %>% group_by(readout, group, mouse_id, first_admin, second_admin, plot_label) %>%
    group_modify(function(d, key) {
      d <- d[order(d$time_h), ]
      missing <- setdiff(expected_times, d$time_h)
      complete <- length(missing) == 0L && nrow(d) == length(expected_times)
      peak <- max(d$value)
      data.frame(n_observed = nrow(d), complete_profile = complete,
                 missing_times_h = paste(missing, collapse = ";"),
                 auc_0_24_ug_h_ml = if (complete) sum(diff(d$time_h) * (head(d$value, -1) + tail(d$value, -1)) / 2) else NA_real_,
                 observed_cmax_ug_ml = peak,
                 first_observed_tmax_h = if (peak > 0) min(d$time_h[d$value == peak]) else NA_real_)
    }) %>% ungroup()
  mean_or_na <- function(x) if (length(x)) mean(x) else NA_real_
  sd_or_na <- function(x) if (length(x) > 1L) sd(x) else NA_real_
  median_or_na <- function(x) if (length(x)) median(x) else NA_real_
  min_or_na <- function(x) if (length(x)) min(x) else NA_real_
  max_or_na <- function(x) if (length(x)) max(x) else NA_real_
  pk_summary <- subject_pk %>% group_by(readout, group, first_admin, second_admin, plot_label) %>%
    group_modify(function(d, key) {
      complete <- d[d$complete_profile, ]
      tmax <- complete$first_observed_tmax_h[!is.na(complete$first_observed_tmax_h)]
      data.frame(n_profiles = nrow(d), n_complete_profiles = nrow(complete), n_incomplete_profiles = sum(!d$complete_profile),
                 auc_0_24_mean = mean_or_na(complete$auc_0_24_ug_h_ml),
                 auc_0_24_sd = sd_or_na(complete$auc_0_24_ug_h_ml),
                 auc_0_24_sem = sd_or_na(complete$auc_0_24_ug_h_ml) / sqrt(nrow(complete)),
                 observed_cmax_mean = mean_or_na(complete$observed_cmax_ug_ml),
                 observed_cmax_sd = sd_or_na(complete$observed_cmax_ug_ml),
                 observed_cmax_sem = sd_or_na(complete$observed_cmax_ug_ml) / sqrt(nrow(complete)),
                 n_observed_tmax = length(tmax), observed_tmax_median_h = median_or_na(tmax),
                 observed_tmax_min_h = min_or_na(tmax), observed_tmax_max_h = max_or_na(tmax))
    }) %>% ungroup() %>% rename(treatment = plot_label)
  methods <- data.frame(
    quantity = c("Concentration", "AUC0-24", "Observed Cmax", "Observed Tmax", "Inference"),
    calculation = c("Mean, SD and SEM across available mice at each readout/group/time; units: ug/mL.",
      "Per-mouse linear trapezoids across 0, 0.5, 4, 8 and 24 h; units: ug*h/mL. Incomplete profiles have missing AUC; missing observations are not imputed and integration is not extended beyond 0-24 h.",
      "Largest measured concentration within 0-24 h; units: ug/mL. Group summaries use complete profiles only.",
      "Earliest sampled time with the observed maximum; units: h. Undefined when all concentrations are zero. Group summary: median and range among complete profiles with a positive observed maximum.",
      "Descriptive analysis. No repeated-measures model, pairwise hypothesis test, fitted half-life or multiplicity adjustment."))
  write.csv(summary_data, file.path(output_dir, "Fig_04E-F_summary.csv"), row.names = FALSE)
  write.csv(concentration_summary, file.path(output_dir, "Fig_04E-F_concentration_statistics.csv"), row.names = FALSE)
  write.csv(subject_pk, file.path(output_dir, "Fig_04E-F_individual_PK.csv"), row.names = FALSE)
  write.csv(pk_summary, file.path(output_dir, "Fig_04E-F_PK_statistics.csv"), row.names = FALSE)
  write.csv(methods, file.path(output_dir, "Fig_04E-F_statistics_methods.csv"), row.names = FALSE)
  list(tables = list(group_summary = summary_data, measurements = plot_data,
                    concentration_statistics = concentration_summary, individual_PK = subject_pk,
                    PK_statistics = pk_summary, statistics_methods = methods),
       sizes = list("Fig. 4e-f - sequential serum profiles" = c(11, 5.2)))
}

print_results <- function(result) {
  fmt <- function(x) vapply(x, function(value) if (is.finite(value))
    format(signif(value, 3), trim = TRUE, scientific = FALSE) else "NA", character(1))
  pair <- function(mean, spread) paste0(fmt(mean), " +/- ", fmt(spread))
  cat("\nFig. 4e,f | Sequential-administration serum concentrations\n")
  cat("Concentrations: mean +/- SEM in ug/mL; parentheses show observed n.\n")
  for (analyte in c("SLC-3010", "hnTCB2")) {
    x <- result$tables$concentration_statistics
    x <- x[as.character(x$readout) == analyte, ]
    start <- x[x$time_h == 0.5, ]
    end <- x[x$time_h == 24, ]
    end <- end[match(as.character(start$group), as.character(end$group)), ]
    tab <- data.frame(Treatment = start$treatment,
                      `0.5 h` = paste0(pair(start$mean, start$sem), " (n=", start$n, ")"),
                      `24 h` = paste0(pair(end$mean, end$sem), " (n=", end$n, ")"), check.names = FALSE)
    cat("\n", analyte, "\n", sep = "")
    print(tab, row.names = FALSE, right = FALSE)
    k <- result$tables$PK_statistics
    k <- k[as.character(k$readout) == analyte, ]
    tmax <- ifelse(k$n_observed_tmax > 0L,
                   paste0(fmt(k$observed_tmax_median_h), " (", fmt(k$observed_tmax_min_h), "-", fmt(k$observed_tmax_max_h), ")"),
                   "not defined")
    tab <- data.frame(Treatment = k$treatment, n = k$n_complete_profiles,
                      `AUC0-24` = pair(k$auc_0_24_mean, k$auc_0_24_sd),
                      Cmax = pair(k$observed_cmax_mean, k$observed_cmax_sd),
                      `Tmax (h)` = tmax, check.names = FALSE)
    cat("Complete profiles: AUC0-24 (ug*h/mL) and observed Cmax (ug/mL), mean +/- SD.\n")
    cat("Observed Tmax: median (range), h; n is the number of complete profiles.\n")
    print(tab, row.names = FALSE, right = FALSE)
  }
  n_missing <- sum(!result$tables$individual_PK$complete_profile)
  cat("\nIncomplete profiles excluded from PK group summaries: ", n_missing,
      "; no missing values were imputed.\n", sep = "")
  cat("All time points, individual PK values and calculation definitions are saved in the CSV files.\n")
  cat("Results: ", result$output_dir, "\n", sep = "")
  invisible(result)
}

plot_results <- function(result) {
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(patchwork))
  group_colors <- c(G1 = "#8A8A8A", G2 = "#2F55FF", G3 = "#E84A5F", G4 = "#008A3D", G5 = "#7A1E76")
  group_levels <- names(group_colors)
  treatments <- unique(result$tables$measurements[c("group", "plot_label")])
  group_labels <- setNames(treatments$plot_label, as.character(treatments$group))
  make_panel <- function(analyte, title) {
    d <- result$tables$group_summary
    d <- d[as.character(d$readout) == analyte, ]
    ggplot(d, aes(time_h, mean, colour = group, group = group)) +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.42, linewidth = 0.45, alpha = 0.92) +
      geom_line(linewidth = 0.65, lineend = "round") + geom_point(size = 2.2, stroke = 0) +
      scale_colour_manual(values = group_colors, breaks = group_levels, labels = group_labels, drop = FALSE) +
      scale_x_continuous(breaks = seq(0, 25, by = 5), expand = expansion(mult = c(0.012, 0.018))) +
      scale_y_continuous(breaks = seq(0, 40, by = 10), expand = expansion(mult = c(0, 0.045))) +
      coord_cartesian(xlim = c(0, 25), ylim = c(0, 40), clip = "off") +
      labs(title = title, x = "Hours after second administration",
           y = paste0("Serum ", analyte, " (\u00b5g/mL)"),
           colour = "First administration \u2192 second administration") +
      theme_classic(base_size = 11, base_family = "Arial") +
      theme(axis.title = element_text(size = 11), axis.text = element_text(size = 10, colour = "black"),
            axis.line = element_line(linewidth = 0.45), axis.ticks = element_line(linewidth = 0.45),
            plot.title = element_text(size = 13, face = "bold", margin = margin(b = 8)),
            plot.margin = margin(7, 12, 7, 7), legend.position = "bottom",
            legend.title = element_text(size = 10), legend.text = element_text(size = 10),
            legend.key.width = grid::unit(18, "pt"), legend.key.height = grid::unit(10, "pt"),
            legend.margin = margin(t = 5), legend.box.spacing = grid::unit(5, "pt")) +
      guides(colour = guide_legend(nrow = 2, byrow = TRUE, title.position = "top", title.hjust = 0.5))
  }
  combined <- (make_panel("SLC-3010", "e  SLC-3010") + make_panel("hnTCB2", "f  hnTCB2") +
    plot_layout(nrow = 1, widths = c(1, 1), guides = "collect")) &
    theme(legend.position = "bottom", legend.justification = "center", legend.box = "vertical")
  list("Fig. 4e-f - sequential serum profiles" = combined)
}

# Run the analysis when the file is sourced.
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
    paths <- c(frames[basename(frames) == "Fig_04E-F_sequential_serum_PK.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_04E-F_sequential_serum_PK.R"], join(getwd(), "Fig_04E-F_sequential_serum_PK.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_04E-F_sequential_serum_PK.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply at most one workbook or directory path.", call. = FALSE)
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
      run <- join(base, "Fig_04E-F_sequential_serum_PK", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
    result <- analyze_data(path, out)
    c(result, list(input = list(path = path, sha256 = sha, status = status), output_dir = out, tables_only = tables_only))
  })
  report_lines <- capture.output(print_results(analysis_result))
  cat(paste(report_lines, collapse = "\n"), "\n", sep = "")
  writeLines(report_lines, file.path(analysis_result$output_dir, "analysis_summary.txt"), useBytes = TRUE)
  if (!analysis_result$tables_only) {
    analysis_result$plots <- tryCatch(local({
      needed <- c("dplyr", "ggplot2", "patchwork")
      missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
      if (length(missing)) stop("Install these plotting packages: ", paste(missing, collapse = ", "))
      if (interactive() && identical(Sys.getenv("RSTUDIO"), "1") && !identical(names(grDevices::dev.cur()), "RStudioGD"))
        grDevices::dev.new(noRStudioGD = FALSE)
      measure <- !interactive() && grDevices::dev.cur() == 1L
      if (measure) grDevices::pdf(NULL)
      plots <- tryCatch(plot_results(analysis_result), finally = if (measure) grDevices::dev.off())
      for (i in seq_along(plots)) {
        size <- analysis_result$sizes[[names(plots)[i]]]
        if (is.null(size)) size <- c(11, 5.2)
        p <- plots[[i]]
        if (is.function(p)) p <- grid::grid.grabExpr(p(), width = size[1L], height = size[2L])
        if (interactive() && grid::is.grob(p)) { grid::grid.newpage(); grid::grid.draw(p) } else if (interactive()) print(p)
        ggplot2::ggsave(file.path(analysis_result$output_dir, sprintf("plot_%02d.png", i)), p,
                        width = size[1L], height = size[2L], units = "in", dpi = 220, bg = "white")
      }
      plots
    }), error = function(e) stop("Analysis tables saved in ", analysis_result$output_dir, ". Plotting failed: ", conditionMessage(e), call. = FALSE))
    figure_result <- analysis_result
  }
}
