# Fig. 2a and Supplementary Fig. 3a | Transferred and host CD8 counts
# Summarize the individual mouse counts for each cell population and treatment and retain
# sample sizes and individual values in millions of cells.
#
# Input: put Source_Data.xlsx beside this script or in the repository's data/ folder.
# Numerical packages: digest, dplyr, readxl, stringr, tibble.
# Plotting packages: dplyr, ggplot2, stringr.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$group_summary
#   analysis_result$output_dir   # Folder containing Fig_02A_Suppl_Fig_03A_summary.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Figure 2a and Supplementary Figure 3a: transferred OT-1 and host CD8 T cell counts.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  library(dplyr)
  library(stringr)

  treatment_levels <- c("Vehicle", "OVA", "BetaIL-2(16.6ug)", "SLC-3010(18ug)", "OVA+BetaIL-2(16.6ug)", "OVA+BetaIL-2(33.2ug)", 
      "OVA+SLC-3010(18ug)", "OVA+SLC-3010(36ug)")

  group_labels <- paste0("G", seq_along(treatment_levels))

  treatment_display_labels <- c("Vehicle", "OVA", "βIL-2 (16.6 µg)", "SLC-3010 (18 µg)", "OVA + βIL-2 (16.6 µg)", 
      "OVA + βIL-2 (33.2 µg)", "OVA + SLC-3010 (18 µg)", "OVA + SLC-3010 (36 µg)")

  stopifnot(length(treatment_levels) == length(treatment_display_labels))

  treatment_axis <- tibble::tibble(treatment = treatment_levels, group_label = group_labels, treatment_label = treatment_display_labels, 
      treatment_order = seq_along(treatment_levels))

  readout_levels <- c("Thy1.1+CD8", "Thy1.1-CD8")

  drug_family_levels <- c("Vehicle", "OVA", "Beta IL-2", "SLC-3010")

  transferred <- readxl::read_excel(workbook, sheet = "Figure 2", range = "A7:D29", .name_repair = "minimal")
  host <- readxl::read_excel(workbook, sheet = "Supplementary Figure 3", range = "A7:D29", .name_repair = "minimal")
  raw_data <- dplyr::bind_rows(
      dplyr::mutate(transferred, readout = "Thy1.1+CD8"),
      dplyr::mutate(host, readout = "Thy1.1-CD8"))
  if (any(!raw_data$group %in% group_labels)) stop("Unknown group in the CD8 count data.")
  # The source group codes use the same treatment-axis order as the figure.
  raw_data <- raw_data %>% dplyr::transmute(mouse_id = mouse_id,
      treatment = treatment_levels[match(group, group_labels)], readout = readout,
      value = cell_count_1e6)

  required_cols <- c("mouse_id", "treatment", "readout", "value")

  missing_cols <- setdiff(required_cols, colnames(raw_data))

  if (length(missing_cols) > 0) {
      stop("Input file is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }

  plot_data <- raw_data %>% mutate(mouse_id = as.character(mouse_id), readout = as.character(readout),
      value = as.numeric(value)) %>% filter(!is.na(treatment), !is.na(readout), !is.na(value)) %>% left_join(treatment_axis, 
      by = "treatment") %>% mutate(readout = factor(readout, levels = readout_levels), drug_family = case_when(treatment == 
      "Vehicle" ~ "Vehicle", treatment == "OVA" ~ "OVA", str_detect(treatment, "BetaIL-2") ~ "Beta IL-2", str_detect(treatment, 
      "SLC-3010") ~ "SLC-3010", TRUE ~ NA_character_), drug_family = factor(drug_family, levels = drug_family_levels))

  if (any(is.na(plot_data$treatment_order))) {
      unknown_treatments <- unique(plot_data$treatment[is.na(plot_data$treatment_order)])
      stop("Unknown treatment label(s) found in input data: ", paste(unknown_treatments, collapse = ", "))
  }

  if (any(is.na(plot_data$readout))) {
      unknown_readouts <- unique(raw_data$readout[!raw_data$readout %in% readout_levels])
      stop("Unknown readout label(s) found in input data: ", paste(unknown_readouts, collapse = ", "))
  }

  plot_data <- plot_data %>% arrange(readout, treatment_order, suppressWarnings(as.numeric(mouse_id)), mouse_id)

  summary_data <- plot_data %>% group_by(readout, treatment, group_label, treatment_label, treatment_order) %>% summarise(n = n(), 
      mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE), sem = sd/sqrt(n), median = median(value, na.rm = TRUE), 
      q1 = quantile(value, 0.25, na.rm = TRUE, names = FALSE), q3 = quantile(value, 0.75, na.rm = TRUE, names = FALSE), 
      min = min(value, na.rm = TRUE), max = max(value, na.rm = TRUE), .groups = "drop") %>% mutate(sd = if_else(is.na(sd), 
      0, sd), sem = if_else(is.na(sem), 0, sem))

  write.csv(summary_data, file = file.path(output_dir, "Fig_02A_Suppl_Fig_03A_summary.csv"), row.names = FALSE)

  individuals <- plot_data %>% select(mouse_id, readout, treatment, group_label, treatment_label, value)
  utils::write.csv(individuals, file.path(output_dir, "Fig_02A_Suppl_Fig_03A_individual_counts.csv"), row.names = FALSE)
  list(tables = list(group_summary = summary_data, individual_counts = individuals),
       measurements = plot_data, treatment_key = treatment_axis)
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result) {
  library(dplyr)
  library(ggplot2)
  library(stringr)

  base_family <- "Arial"
  font_base_pt <- 6.2
  font_panel_title_pt <- 7
  font_axis_title_pt <- 6.3
  font_x_tick_pt <- 6.2
  font_y_tick_pt <- 6
  show_panel_title <- TRUE
  add_group_separators <- FALSE
  show_group_legend <- TRUE
  legend_title <- NULL

  legend_position <- c(0.015, 0.99)

  legend_text_size <- 6
  legend_key_height_pt <- 6.3
  legend_key_width_pt <- 6
  legend_point_size <- 1.25
  transferred_ymax <- 25
  transferred_ybreak_step <- 5
  host_ymax <- 125
  host_ybreak_step <- 25

  treatment_axis <- result$treatment_key
  treatment_levels <- treatment_axis$treatment
  group_labels <- treatment_axis$group_label
  treatment_display_labels <- treatment_axis$treatment_label
  group_legend_labels <- paste0(group_labels, ": ", treatment_display_labels)
  plot_data <- result$measurements

  drug_family_levels <- c("Vehicle", "OVA", "Beta IL-2", "SLC-3010")

  drug_colors <- c(Vehicle = "#4D4D4D", OVA = "#000000", `Beta IL-2` = "#3B5FD9", `SLC-3010` = "#E53935")

  drug_fills <- c(Vehicle = "white", OVA = "#000000", `Beta IL-2` = "#3B5FD9", `SLC-3010` = "#E53935")

  drug_shapes <- c(Vehicle = 21, OVA = 16, `Beta IL-2` = 17, `SLC-3010` = 15)

  treatment_drug_family <- dplyr::case_when(treatment_levels == "Vehicle" ~ "Vehicle", treatment_levels == "OVA" ~ 
      "OVA", stringr::str_detect(treatment_levels, "BetaIL-2") ~ "Beta IL-2", stringr::str_detect(treatment_levels, 
      "SLC-3010") ~ "SLC-3010", TRUE ~ NA_character_)

  if (anyNA(treatment_drug_family)) {
      stop("Each treatment must have a drug-family color and shape mapping.")
  }

  treatment_shapes <- stats::setNames(unname(drug_shapes[treatment_drug_family]), treatment_levels)

  treatment_colors <- unname(drug_colors[treatment_drug_family])

  treatment_fills <- unname(drug_fills[treatment_drug_family])

  boxplot_width <- 0.44
  boxplot_linewidth <- 0.28
  boxplot_alpha <- 0.16
  point_jitter_width <- 0.1
  point_size <- 1.3
  point_stroke <- 0.34
  transferred_y_title <- "Thy1.1⁺CD8 (×10⁶)"
  host_y_title <- "Thy1.1⁻CD8 (×10⁶)"
  transferred_panel_title <- "Thy1.1⁺CD8⁺"
  host_panel_title <- "Thy1.1⁻CD8⁺"

  plot_data <- plot_data %>% arrange(readout, treatment_order, suppressWarnings(as.numeric(mouse_id)), mouse_id) %>% 
      group_by(readout, treatment) %>% mutate(replicate_index = row_number(), n_in_group = n(), jitter_x = if_else(n_in_group <= 
      1, 0, (replicate_index - (n_in_group + 1)/2) * (2 * point_jitter_width/pmax(n_in_group - 1, 1))), x_point = treatment_order + 
      jitter_x) %>% ungroup()

  theme_inside_legend <- function() {
      if (!isTRUE(show_group_legend)) 
          return(theme(legend.position = "none"))
      if (utils::packageVersion("ggplot2") >= "3.5.0") 
          theme(legend.position = "inside", legend.position.inside = legend_position, legend.justification = c(0, 
              1), legend.justification.inside = c(0, 1), legend.key.spacing.y = grid::unit(0, "pt"))
      else theme(legend.position = legend_position, legend.justification = c(0, 1))
  }

  count_theme <- function(base_size = font_base_pt, base_family = "Arial") theme_classic(base_size = base_size, base_family = base_family) + 
      theme(plot.title = element_text(size = font_panel_title_pt, face = "bold", color = "black", hjust = 0.5, margin = margin(b = 0.8)), 
          axis.title.x = element_blank(), axis.title.y = element_text(size = font_axis_title_pt, color = "black", 
              angle = 90, hjust = 0.5, vjust = 0.5, margin = margin(r = 2)), axis.text.x = element_text(size = font_x_tick_pt, 
              color = "black", angle = 0, hjust = 0.5, vjust = 1, margin = margin(t = 1.5)), axis.text.y = element_text(size = font_y_tick_pt, 
              color = "black", margin = margin(r = 0.7)), axis.line = element_line(color = "black", linewidth = 0.25), 
          axis.ticks = element_line(color = "black", linewidth = 0.25), axis.ticks.length = unit(1.1, "pt"), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), legend.direction = "vertical", legend.title = element_blank(), legend.text = element_text(size = legend_text_size, 
              color = "black", hjust = 0, margin = margin(l = 1.5)), legend.key = element_blank(), legend.key.height = grid::unit(legend_key_height_pt, 
              "pt"), legend.key.width = grid::unit(legend_key_width_pt, "pt"), legend.background = element_blank(), 
          legend.box.background = element_blank(), legend.margin = margin(0, 0.5, 0, 0.5, unit = "pt"), legend.box.margin = margin(0, 
              0, 0, 0, unit = "pt"), legend.spacing.y = grid::unit(0, "pt"), plot.margin = margin(t = 0.6, r = 0.8, 
              b = 0.8, l = 0.8, unit = "pt")) + theme_inside_legend()

  make_readout_plot <- function(data, readout_value, y_limits, y_breaks, y_title, panel_title, base_size = font_base_pt) {
      panel_data <- data %>% filter(as.character(readout) == readout_value)
      if (nrow(panel_data) == 0L) 
          stop("No data found for readout: ", readout_value)
      if (any(!is.finite(panel_data$value))) stop("Non-finite counts for ", readout_value)
      if (any(panel_data$value < y_limits[1] | panel_data$value > y_limits[2])) {
          y_limits <- range(c(y_limits, panel_data$value))
          y_breaks <- pretty(y_limits)
      }
      p <- ggplot(panel_data) + (if (isTRUE(add_group_separators)) 
          geom_vline(xintercept = c(2.5, 4.5), linewidth = 0.24, color = "black", alpha = 0.14)) + geom_boxplot(aes(x = treatment_order, 
          y = value, group = treatment, fill = drug_family), width = boxplot_width, linewidth = boxplot_linewidth, 
          color = "black", alpha = boxplot_alpha, outlier.shape = NA, show.legend = FALSE) + geom_point(aes(x = x_point, 
          y = value, color = drug_family, fill = drug_family, shape = treatment), size = point_size, stroke = point_stroke, 
          alpha = 1, show.legend = c(colour = FALSE, fill = FALSE, shape = TRUE)) + scale_x_continuous(breaks = treatment_axis$treatment_order, 
          labels = treatment_axis$group_label, limits = c(0.55, length(treatment_levels) + 0.55), expand = expansion(mult = c(0, 
              0))) + scale_y_continuous(limits = y_limits, breaks = y_breaks, expand = expansion(mult = c(0, 0.035))) + 
          scale_color_manual(values = drug_colors, breaks = drug_family_levels, guide = "none", drop = FALSE) + scale_fill_manual(values = drug_fills, 
          breaks = drug_family_levels, guide = "none", drop = FALSE) + scale_shape_manual(name = legend_title, values = treatment_shapes, 
          limits = treatment_levels, breaks = treatment_levels, labels = group_legend_labels, drop = FALSE, guide = guide_legend(ncol = 1, 
              override.aes = list(colour = treatment_colors, fill = treatment_fills, shape = unname(treatment_shapes), 
                  size = legend_point_size, stroke = point_stroke, alpha = 1))) + labs(x = NULL, y = y_title, title = if (isTRUE(show_panel_title)) 
          panel_title
      else NULL) + coord_cartesian(clip = "off") + count_theme(base_size = base_size, base_family = base_family)
      p
  }

  transferred_plot <- make_readout_plot(data = plot_data, readout_value = "Thy1.1+CD8", y_limits = c(0, transferred_ymax), y_breaks = seq(0, 
      transferred_ymax, by = transferred_ybreak_step), y_title = transferred_y_title, panel_title = transferred_panel_title, base_size = font_base_pt)

  host_plot <- make_readout_plot(data = plot_data, readout_value = "Thy1.1-CD8", y_limits = c(0, host_ymax), y_breaks = seq(0, 
      host_ymax, by = host_ybreak_step), y_title = host_y_title, panel_title = host_panel_title, base_size = font_base_pt)

  list("Fig. 2a - transferred OT-1 CD8 T cells" = transferred_plot,
                    "Supplementary Fig. 3a - host CD8 T cells" = host_plot)
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
  cat("Fig. 2a and Supplementary Fig. 3a | Transferred and host CD8 counts\nDescriptive statistics across individual mice; no additional hypothesis tests.\n")
  show("group_summary", "Counts by treatment (millions of cells)", c("readout", "treatment_label", "n", "mean", "sd", "sem", "median", "q1", "q3"))
  cat("\nComplete result tables: ", paste(names(result$tables), collapse = ", "), "\n", sep = "")
  cat("Results directory: ", result$output_dir, "\n", sep = "")
  invisible(result)
}

# Source calculates the results, prints the analysis summary and displays the plots.
if (isTRUE(getOption("SLC3010.autorun", TRUE))) {
  analysis_result <- figure_result <- NULL
  analysis_result <- local({
    needed <- c("digest", "dplyr", "readxl", "stringr", "tibble")
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
    paths <- c(frames[basename(frames) == "Fig_02A_Suppl_Fig_03A_OT1_Host_CD8.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_02A_Suppl_Fig_03A_OT1_Host_CD8.R"], join(getwd(), "Fig_02A_Suppl_Fig_03A_OT1_Host_CD8.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_02A_Suppl_Fig_03A_OT1_Host_CD8.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply one workbook path or directory, or use Source in RStudio.", call. = FALSE)
    path <- path.expand(if (length(args)) args[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR"))
    book <- "Source_Data.xlsx"
    sheets <- c("Figure 2", "Supplementary Figure 3")
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
      run <- join(base, "Fig_02A_Suppl_Fig_03A_OT1_Host_CD8", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
      needed <- c("dplyr", "ggplot2", "stringr")
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
        if (is.null(size)) size <- c(9, 5)
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
