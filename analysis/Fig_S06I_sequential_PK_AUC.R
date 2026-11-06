# Supplementary Fig. 6i | Individual-mouse sequential serum exposure
# Integrate observed concentration-time profiles with the linear trapezoidal rule per mouse
# and analyte, then summarize those AUCs by treatment.
# Check each mouse's first and last measured times and the number of measurements before
# comparing AUCs. The observation windows can differ.
#
# Input: put Source_Data.xlsx beside this script or in the repository's data/ folder.
# Numerical packages: digest, dplyr, readxl.
# Plotting packages: dplyr, ggplot2, scales.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$mouse_auc
#   analysis_result$output_dir   # Folder containing Fig_S06I_mouse_AUC.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Supplementary Fig. 6i: Sequential serum profiles and individual-mouse trapezoidal AUC summaries.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  library(dplyr)

  group_levels <- c("PBS->PBS", "hnTCB2->IL-2", "SLC-3010->IL-2", "PBS->hnTCB2", "PBS->SLC-3010")

  treatment_id_labels <- c(`PBS->PBS` = "G1", `hnTCB2->IL-2` = "G2", `SLC-3010->IL-2` = "G3", `PBS->hnTCB2` = "G4", `PBS->SLC-3010` = "G5")

  readout_levels <- c("SLC-3010", "hnTCB2")

  readout_labels <- c(`SLC-3010` = "SLC-3010", hnTCB2 = "hnTCB2")

  PK_raw <- readxl::read_excel(workbook, sheet = "Supplementary Figure 6", range = "A301:E601", .name_repair = "minimal")

  names(PK_raw) <- trimws(names(PK_raw))

  required_cols <- c("Group", "ID", "Detection", "Time", "Conc")

  missing_cols <- setdiff(required_cols, names(PK_raw))

  if (length(missing_cols) > 0) {
      stop("Input file is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }

  normalize_group <- function(x) {
      x <- as.character(x)
      x <- gsub("\\s*\\u2192\\s*", "->", x)
      x <- gsub("\\s*->\\s*", "->", x)
      x <- gsub("hIL-2", "IL-2", x, fixed = TRUE)
      x
  }

  PK <- PK_raw %>% dplyr::select(dplyr::all_of(required_cols)) %>% dplyr::mutate(Group = normalize_group(Group), 
      ID = as.factor(ID), Detection = as.character(Detection), Time = as.numeric(Time), Conc = as.numeric(Conc)) %>% 
      dplyr::filter(!is.na(Group), !is.na(ID), !is.na(Detection), !is.na(Time), !is.na(Conc)) %>% dplyr::filter(Group %in% 
      group_levels) %>% dplyr::mutate(Group = factor(Group, levels = group_levels), Treatment_ID = factor(as.character(Group), 
      levels = group_levels, labels = treatment_id_labels), Detection = factor(Detection, levels = readout_levels, 
      labels = readout_labels)) %>% dplyr::filter(!is.na(Detection))

  if (nrow(PK) == 0) {
      stop("No data remained after filtering. Check Group and Detection in the Supplementary Figure 6i Source Data block.")
  }

  time_levels_numeric <- c(0, 0.5, 4, 8, 24)

  observed_times <- sort(unique(PK$Time))

  if (!all(observed_times %in% time_levels_numeric)) {
      time_levels_numeric <- observed_times
  }

  format_time_label <- function(x) {
      out <- ifelse(abs(x - round(x)) < 1e-08, as.character(as.integer(round(x))), as.character(x))
      out
  }

  time_levels_label <- format_time_label(time_levels_numeric)

  PK <- PK %>% dplyr::mutate(Time_factor = factor(format_time_label(Time), levels = time_levels_label))

  PK_summary <- PK %>% dplyr::group_by(Group, Treatment_ID, Detection, Time, Time_factor) %>% dplyr::summarise(n = dplyr::n(), 
      mean_conc = mean(Conc, na.rm = TRUE), sd_conc = stats::sd(Conc, na.rm = TRUE), .groups = "drop") %>% dplyr::mutate(sd_conc = ifelse(is.na(sd_conc), 
      0, sd_conc), ymin = pmax(mean_conc - sd_conc, 0), ymax = mean_conc + sd_conc)

  auc_trapezoid <- function(time, conc) {
      ok <- is.finite(time) & is.finite(conc)
      time <- time[ok]
      conc <- conc[ok]
      if (length(time) < 2) 
          return(NA_real_)
      ord <- order(time)
      time <- time[ord]
      conc <- conc[ord]
      sum(diff(time) * (head(conc, -1) + tail(conc, -1))/2)
  }

  # Integrate between the first and last available observations for each mouse.
  mouse_auc <- PK %>% group_by(Group, Treatment_ID, Detection, ID) %>%
      summarise(n_observations = n(), first_h = min(Time), last_h = max(Time),
                observed_hours = paste(sort(Time), collapse = ", "),
                AUC = auc_trapezoid(Time, Conc), .groups = "drop")
  group_auc <- mouse_auc %>% group_by(Group, Treatment_ID, Detection) %>%
      summarise(n = n(), mean_AUC = mean(AUC, na.rm = TRUE), sd_AUC = sd(AUC, na.rm = TRUE),
                n_calculable = sum(is.finite(AUC)), .groups = "drop") %>%
      mutate(sd_AUC = ifelse(is.na(sd_AUC), 0, sd_AUC))

  write.csv(PK_summary, file.path(output_dir, "Fig_S06I_PK_summary.csv"), row.names = FALSE)

  write.csv(mouse_auc, file.path(output_dir, "Fig_S06I_mouse_AUC.csv"), row.names = FALSE)
  write.csv(group_auc, file.path(output_dir, "Fig_S06I_AUC_summary.csv"), row.names = FALSE)
  list(tables = list(mouse_auc = mouse_auc, group_auc = group_auc,
                     concentration_summary = PK_summary), times = time_levels_label)
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result, display = FALSE) {
  library(ggplot2)
  library(dplyr)

  base_family <- "Arial"

  group_colors <- c(`PBS->PBS` = "#8A8A8A", `hnTCB2->IL-2` = "#2F55FF", `SLC-3010->IL-2` = "#E84A5F", `PBS->hnTCB2` = "#008A3D", 
      `PBS->SLC-3010` = "#7A1E76")

  group_levels <- names(group_colors)

  arrow_symbol <- "→"

  group_labels <- c(`PBS->PBS` = paste0("G1 PBS", arrow_symbol, "PBS"), `hnTCB2->IL-2` = paste0("G2 hnTCB2", arrow_symbol, 
      "IL-2"), `SLC-3010->IL-2` = paste0("G3 SLC-3010", arrow_symbol, "IL-2"), `PBS->hnTCB2` = paste0("G4 PBS", arrow_symbol, 
      "hnTCB2"), `PBS->SLC-3010` = paste0("G5 PBS", arrow_symbol, "SLC-3010"))

  treatment_id_labels <- c(`PBS->PBS` = "G1", `hnTCB2->IL-2` = "G2", `SLC-3010->IL-2` = "G3", `PBS->hnTCB2` = "G4", `PBS->SLC-3010` = "G5")

  readout_levels <- c("SLC-3010", "hnTCB2")

  readout_labels <- c(`SLC-3010` = "SLC-3010", hnTCB2 = "hnTCB2")

  PK_summary <- result$tables$concentration_summary
  time_levels_label <- result$times
  AUC_summary <- result$tables$group_auc %>% mutate(
      auc_label = paste0("AUC = ", scales::number(mean_AUC, accuracy = 0.1)),
      Time_factor = factor(tail(time_levels_label, 1), levels = time_levels_label), y_label = 43)

  theme_suppl_pk <- function(base_size = 6.2, base_family = "Arial") ggplot2::theme_classic(base_size = base_size, 
      base_family = base_family) %+replace% ggplot2::theme(axis.title.x = ggplot2::element_text(color = "black", 
      size = base_size + 0.15, margin = ggplot2::margin(t = 3.2)), axis.title.y = ggplot2::element_text(color = "black", 
      size = base_size + 0.15, angle = 90, hjust = 0.5, vjust = 0.5, margin = ggplot2::margin(r = 4.2)), axis.text.x = ggplot2::element_text(color = "black", 
      size = base_size - 0.45, margin = ggplot2::margin(t = 1.1)), axis.text.y = ggplot2::element_text(color = "black", 
      size = base_size - 0.55, margin = ggplot2::margin(r = 1.1)), axis.line = ggplot2::element_line(color = "black", 
      linewidth = 0.3), axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3), axis.ticks.length = grid::unit(1.15, 
      "pt"), panel.border = ggplot2::element_rect(fill = NA, color = "#AFAFAF", linewidth = 0.32), panel.grid = ggplot2::element_blank(), 
      panel.spacing.x = grid::unit(2, "pt"), panel.spacing.y = grid::unit(2.6, "pt"), strip.background = ggplot2::element_rect(fill = "#F3F3F3", 
          color = "#BDBDBD", linewidth = 0.28), strip.text.x = ggplot2::element_text(color = "black", size = base_size + 
          0.35, face = "bold", margin = ggplot2::margin(t = 3.2, b = 2.4)), strip.text.y.right = ggplot2::element_text(color = "black", 
          size = base_size - 0.45, face = "plain", angle = 0, hjust = 0.5, vjust = 0.5, margin = ggplot2::margin(l = 1, 
              r = 1)), legend.position = "top", legend.justification = c(0, 1), legend.direction = "horizontal", 
      legend.box = "vertical", legend.box.just = "left", legend.title = ggplot2::element_blank(), legend.text = ggplot2::element_text(size = base_size - 
          1.35, color = "black", hjust = 0), legend.key = ggplot2::element_blank(), legend.key.width = grid::unit(7, 
          "pt"), legend.key.height = grid::unit(4.8, "pt"), legend.spacing.x = grid::unit(4, "pt"), legend.spacing.y = grid::unit(1, 
          "pt"), legend.margin = ggplot2::margin(t = 0.5, r = 0, b = 0.5, l = 0, unit = "pt"), legend.box.margin = ggplot2::margin(t = 3, 
          r = 0, b = 1.8, l = 0, unit = "pt"), legend.box.spacing = grid::unit(0, "pt"), plot.margin = ggplot2::margin(t = 10, 
          r = 1, b = 1, l = 3, unit = "pt"), plot.background = ggplot2::element_rect(fill = "white", color = NA), 
      panel.background = ggplot2::element_rect(fill = "white", color = NA))

  auc_label_size_mm <- if (display) 11 / ggplot2::.pt else 1.95

  sequential_pk_auc_plot <- ggplot2::ggplot(PK_summary, ggplot2::aes(x = Time_factor, y = mean_conc, color = Group, fill = Group, 
      group = Group)) + ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = mean_conc), alpha = 0.18, color = NA, 
      show.legend = FALSE) + ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax), width = 0.13, linewidth = 0.3, 
      alpha = 0.95, show.legend = FALSE) + ggplot2::geom_line(linewidth = 0.46, lineend = "round", alpha = 0.98) + 
      ggplot2::geom_point(size = 1.15, stroke = 0, alpha = 0.98) + ggplot2::geom_text(data = AUC_summary, ggplot2::aes(x = Time_factor, 
      y = y_label, label = auc_label), inherit.aes = FALSE, hjust = 1.03, vjust = 1.05, size = auc_label_size_mm, 
      family = base_family, color = "black") + ggplot2::facet_grid(rows = ggplot2::vars(Treatment_ID), cols = ggplot2::vars(Detection)) + 
      ggplot2::scale_color_manual(values = group_colors, breaks = group_levels, labels = group_labels, drop = FALSE) + 
      ggplot2::scale_fill_manual(values = group_colors, breaks = group_levels, labels = group_labels, drop = FALSE) + 
      ggplot2::scale_x_discrete(drop = FALSE) + ggplot2::scale_y_continuous(limits = c(0, 45), breaks = seq(0, 40, 
      by = 10), labels = scales::label_number(accuracy = 1), expand = ggplot2::expansion(mult = c(0, 0.02))) + ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, 
      byrow = TRUE, override.aes = list(linewidth = 0.55, size = 1.35, alpha = 1)), fill = "none") + ggplot2::labs(x = expression("Time (hours) after 2"^"nd" ~ 
      "injection"), y = paste0("Concentration (", "µ", "g/mL)")) + ggplot2::coord_cartesian(ylim = c(0, 45), clip = "off") + 
      theme_suppl_pk(base_size = if (display) 11 else 6.2, base_family = base_family)

  list("Supplementary Fig. 6i - individual PK profiles" = sequential_pk_auc_plot)
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
  cat("Supplementary Fig. 6i | Individual-mouse sequential serum exposure\nLinear trapezoidal integration uses each mouse's first-to-last observed interval; windows may differ. No additional hypothesis tests are applied.\n")
  show("group_auc", "Group summaries of individual-mouse AUC (micrograms-hour/mL)", c("Group", "Detection", "n", "n_calculable", "mean_AUC", "sd_AUC"))
  show("mouse_auc", "Individual-mouse AUC and observed integration windows", c("Group", "Detection", "ID", "n_observations", "first_h", "last_h", "AUC"))
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
    paths <- c(frames[basename(frames) == "Fig_S06I_sequential_PK_AUC.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_S06I_sequential_PK_AUC.R"], join(getwd(), "Fig_S06I_sequential_PK_AUC.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_S06I_sequential_PK_AUC.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply one workbook path or directory, or use Source in RStudio.", call. = FALSE)
    path <- path.expand(if (length(args)) args[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR"))
    book <- "Source_Data.xlsx"
    sheets <- c("Supplementary Figure 6")
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
      run <- join(base, "Fig_S06I_sequential_PK_AUC", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
