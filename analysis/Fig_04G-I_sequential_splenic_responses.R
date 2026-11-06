# Fig. 4g-i | Splenic responses after sequential administration
# Summarize each individual-mouse readout within treatment groups and retain the treatment
# schedule and number of contributing observations.
#
# Input: put Source_Data.xlsx beside this script or in the repository's data/ folder.
# Numerical packages: digest, dplyr, readxl, stringr.
# Plotting packages: dplyr, ggbeeswarm, ggplot2, patchwork, scales.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$summary
#   analysis_result$output_dir   # Folder containing Fig_04G-I_summary.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Fig. 4g-i: Splenic MP CD8 T cells, Tregs and their ratio after sequential treatment.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  suppressPackageStartupMessages(library(dplyr))

  group_levels <- c("G1", "G2", "G3", "G4", "G5")

  readout_levels <- c("MP_CD8_T_cells", "Treg", "MP_CD8_T_Treg_ratio")

  raw_data <- readxl::read_excel(workbook, sheet = "Figure 4", range = "A730:F775", .name_repair = "minimal")
  raw_data$plot_label <- paste0(raw_data$first_admin, "→", raw_data$second_admin)
  raw_data <- raw_data[c("mouse_id", "group", "first_admin", "second_admin", "plot_label", "readout", "value")]

  colnames(raw_data) <- stringr::str_trim(colnames(raw_data))

  required_cols <- c("mouse_id", "group", "first_admin", "second_admin", "plot_label", "readout", "value")

  missing_cols <- setdiff(required_cols, colnames(raw_data))

  if (length(missing_cols) > 0) {
      stop("Input file is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }

  measurements <- raw_data %>% mutate(mouse_id = as.character(mouse_id), group = factor(as.character(group), levels = group_levels), 
      first_admin = as.character(first_admin), second_admin = as.character(second_admin), plot_label = as.character(plot_label), 
      readout = as.character(readout), value = as.numeric(value)) %>% filter(!is.na(group), !is.na(readout), !is.na(value)) %>% 
      filter(readout %in% readout_levels) %>% mutate(readout = factor(readout, levels = readout_levels))

  if (nrow(measurements) == 0) {
      stop("No eligible observations remained. Check readout names in the Figure 4g–i Source Data block.")
  }

  missing_readouts <- setdiff(readout_levels, unique(as.character(measurements$readout)))

  if (length(missing_readouts) > 0) {
      warning("These readout(s) were not found in the data: ", paste(missing_readouts, collapse = ", "))
  }

  summary_data <- measurements %>% group_by(readout, group) %>% summarise(n = n(), mean = mean(value, na.rm = TRUE), 
      sd = sd(value, na.rm = TRUE), sem = sd/sqrt(n), .groups = "drop") %>% mutate(readout = factor(as.character(readout), 
      levels = readout_levels))

  decisions <- unique(measurements[c("group", "first_admin", "second_admin")])
  if (anyDuplicated(decisions$group)) stop("Each group needs one treatment schedule.")
  decisions <- decisions[order(decisions$group), ]

  # Preserve the original display: on logarithmic axes, the bars summarize
  # the log10 measurements; quantitative arithmetic means are reported separately.
  view <- measurements %>% mutate(axis_scale = if_else(as.character(readout) == "MP_CD8_T_Treg_ratio", "linear", "log10"),
                           display_value = value)
  log_rows <- view$axis_scale == "log10" & view$value > 0
  view$display_value[log_rows] <- log10(view$value[log_rows])
  view$display_value[view$axis_scale == "log10" & view$value <= 0] <- NA_real_
  display_summary <- view %>% group_by(readout, group, axis_scale) %>%
    summarise(n_display = sum(is.finite(display_value)), centre = mean(display_value, na.rm = TRUE),
              spread = sd(display_value, na.rm = TRUE) / sqrt(n_display), .groups = "drop") %>%
    mutate(lower = pmax(centre - spread, 1e-06), upper = centre + spread,
           across(c(centre, lower, upper), ~ if_else(axis_scale == "log10", 10 ^ .x, .x)))

  write.csv(summary_data, file.path(output_dir, "Fig_04G-I_summary.csv"), row.names = FALSE)
  write.csv(measurements, file.path(output_dir, "Fig_04G-I_individual_measurements.csv"), row.names = FALSE)
  write.csv(display_summary, file.path(output_dir, "Fig_04G-I_display_summary.csv"), row.names = FALSE)
  list(tables = list(summary = summary_data, observations = measurements, treatment_schedule = decisions,
                     display_summary = display_summary))
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result) {
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(scales))
  suppressPackageStartupMessages(library(grid))

  plot_data <- result$tables$observations
  decisions <- result$tables$treatment_schedule
  display_summary <- result$tables$display_summary
  group_levels <- levels(plot_data$group)
  readout_levels <- levels(plot_data$readout)

  base_family <- "Arial"

  group_colors <- c(G1 = "#8A8A8A", G2 = "#2F55FF", G3 = "#E84A5F", G4 = "#008A3D", G5 = "#7A1E76")

  readout_titles <- c(MP_CD8_T_cells = "MP CD8 T cells", Treg = "Treg", MP_CD8_T_Treg_ratio = "MP CD8 T / Treg ratio"
  )

  readout_y_labels <- list(MP_CD8_T_cells = "# of MP CD8 T cells (x10⁶)", Treg = "# of Tregs (x10⁶)", MP_CD8_T_Treg_ratio = "MP CD8 T / Treg ratio")

  y_axis_config <- list(MP_CD8_T_cells = list(scale_type = "log10", limits = c(1, 100), breaks = c(1, 3, 10, 30, 
  100)), Treg = list(scale_type = "log10", limits = c(1, 30), breaks = c(1, 3, 10, 30)), MP_CD8_T_Treg_ratio = list(scale_type = "linear", 
      limits = c(0, 5), breaks = seq(0, 5, by = 1)))

  pt_size_circle <- 3.65
  pt_stroke <- 0.75
  mean_bar_width <- 0.46
  sem_width <- 0.13
  box_width <- 0.52
  box_alpha <- 0.15
  base_font_size <- 10.5
  x_text_size <- 10.5
  y_text_size <- 9.8
  title_size <- 11.2
  condition_tile_row_label_size <- 9
  condition_tile_value_size <- 2.65
  condition_tile_long_value_size <- 2.35
  condition_tile_row_label_margin_r <- 4.5
  condition_tile_row_label_face <- "plain"
  condition_tile_rel_height <- 0.22
  y_label_axis_gap_pt <- 0
  y_title_overlay_x_panel <- 0.116
  y_title_overlay_y_panel_no_tiles <- 0.515
  y_title_overlay_width_panel <- 0.165
  y_title_overlay_height_panel_no_tiles <- 0.86
  y_title_overlay_height_panel_with_tiles <- 0.7

  splenic_response_theme <- function(base_size = base_font_size, base_family_this = base_family, x_text_angle = 0, 
      x_text_hjust = 0.5, x_text_vjust = 0.7, top_margin = 8, bottom_margin = 10, hide_y_title = FALSE) theme_classic(base_size = base_size, 
      base_family = base_family_this) + theme(plot.title = element_text(face = "bold", size = title_size, hjust = 0.5, 
      margin = margin(b = 4)), axis.title.x = element_blank(), axis.title.y = element_text(color = if (isTRUE(hide_y_title)) "transparent" else "black", 
      size = y_text_size, margin = margin(r = y_label_axis_gap_pt)), axis.text.y = element_text(color = "black", 
      size = y_text_size), axis.text.x = element_text(color = "black", size = x_text_size, angle = x_text_angle, 
      hjust = x_text_hjust, vjust = x_text_vjust, margin = margin(t = 3)), axis.line = element_line(color = "black", 
      linewidth = 0.4), axis.ticks = element_line(color = "black", linewidth = 0.4), axis.ticks.length = unit(1.8, 
      "mm"), legend.position = "none", panel.grid = element_blank(), plot.margin = margin(top_margin, 12, bottom_margin, 
      4))

  get_y_scale <- function(readout_name) {
      marks <- display_summary[as.character(display_summary$readout) == readout_name, , drop = FALSE]
      y_cfg <- y_axis_config[[readout_name]]
      if (y_cfg$scale_type == "log10") 
          scale_y_log10(breaks = y_cfg$breaks, labels = label_number(accuracy = 1), expand = expansion(mult = c(0.02, 
              0.08)))
      else scale_y_continuous(breaks = y_cfg$breaks, labels = label_number(accuracy = 1), expand = expansion(mult = c(0.02, 
          0.08)))
  }

  make_splenic_panel <- function(readout_name, show_title = TRUE, x_text_angle = 0, x_text_hjust = 0.5, x_text_vjust = 0.7, 
      hide_y_title = FALSE) {
      df_one <- plot_data %>% filter(as.character(readout) == readout_name) %>% mutate(group = factor(as.character(group), 
          levels = group_levels))
      if (nrow(df_one) == 0) 
          stop("No data for readout: ", readout_name)
      n_by_group <- df_one %>% count(group, name = "n_group")
      df_box <- df_one %>% left_join(n_by_group, by = "group") %>% filter(n_group >= 5)
      marks <- display_summary[as.character(display_summary$readout) == readout_name, , drop = FALSE]
      y_cfg <- y_axis_config[[readout_name]]
      p <- ggplot(df_one, aes(x = group, y = value)) + geom_boxplot(data = df_box, aes(fill = group, color = group), 
          width = box_width, linewidth = 0.38, alpha = box_alpha, outlier.shape = NA, show.legend = FALSE) + geom_errorbar(data = marks, aes(x = group, ymin = lower, ymax = upper),
          inherit.aes = FALSE, width = sem_width, linewidth = 0.45, colour = "black") +
          geom_errorbar(data = marks, aes(x = group, ymin = centre, ymax = centre),
          inherit.aes = FALSE, width = mean_bar_width, linewidth = 0.45, colour = "black") + 
          ggbeeswarm::geom_quasirandom(aes(fill = group, color = group), shape = 21, width = 0.18, size = pt_size_circle, 
              stroke = pt_stroke, alpha = 0.98, varwidth = FALSE, show.legend = FALSE) + scale_fill_manual(values = group_colors, 
          drop = FALSE) + scale_color_manual(values = group_colors, drop = FALSE) + scale_x_discrete(limits = group_levels, 
          drop = FALSE, expand = expansion(add = c(0.6, 0.6))) + get_y_scale(readout_name) + labs(title = if (isTRUE(show_title)) 
          readout_titles[[readout_name]]
      else NA_character_, y = readout_y_labels[[readout_name]]) + splenic_response_theme(x_text_angle = x_text_angle, 
          x_text_hjust = x_text_hjust, x_text_vjust = x_text_vjust, hide_y_title = hide_y_title) + coord_cartesian(ylim = y_cfg$limits, 
          clip = "off")
      if (readout_name %in% c("MP_CD8_T_cells", "Treg")) 
          p <- p + annotation_logticks(sides = "l", outside = FALSE, short = unit(0.6, "mm"), mid = unit(0.95, "mm"), 
              long = unit(1.35, "mm"), linewidth = 0.27)
      p
  }

  condition_table_design <- data.frame(
    condition = rep(c("1st admin.", "2nd admin."), each = nrow(decisions)),
    group = rep(decisions$group, times = 2),
    value = c(decisions$first_admin, decisions$second_admin)
  ) %>% mutate(condition = factor(condition, levels = c("2nd admin.", "1st admin.")),
               fill_hex = case_when(value == "PBS" ~ "#F7F7F7", value == "hnTCB2" ~ "#E8EEFF",
                                    value == "hIL-2" ~ "#F8D8E8", value == "SLC-3010" ~ "#E6E0FF", TRUE ~ "#EFEFEF"),
               text_hex = "#111111",
               value_size_mm = if_else(value == "SLC-3010", condition_tile_long_value_size, condition_tile_value_size))

  make_condition_tile_table <- function() ggplot(condition_table_design, aes(x = group, y = condition)) + geom_tile(aes(fill = fill_hex), 
      width = 0.88, height = 0.78, color = "white", linewidth = 0.35, show.legend = FALSE) + geom_text(aes(label = value, 
      color = text_hex, size = value_size_mm), family = base_family, fontface = "bold", show.legend = FALSE) + scale_fill_identity(guide = "none") + 
      scale_color_identity(guide = "none") + scale_size_identity(guide = "none") + scale_x_discrete(limits = group_levels, 
      drop = FALSE, expand = expansion(add = c(0.6, 0.6)), position = "top") + scale_y_discrete(drop = FALSE) + coord_cartesian(clip = "off") + 
      labs(x = NULL, y = NULL) + theme_void(base_family = base_family) + theme(axis.text.x = element_text(size = 9), axis.text.y = element_text(color = "black", 
      size = condition_tile_row_label_size, face = condition_tile_row_label_face, hjust = 1, margin = margin(r = condition_tile_row_label_margin_r), 
      lineheight = 0.92), axis.ticks = element_blank(), panel.grid = element_blank(), panel.background = element_rect(fill = "white", 
      color = NA), plot.background = element_rect(fill = "white", color = NA), plot.margin = margin(0, 12, 4, 12))

  charts <- patchwork::wrap_plots(lapply(readout_levels, make_splenic_panel), nrow = 1)
  splenic_response_composite <- patchwork::wrap_plots(charts, make_condition_tile_table(),
                                                      ncol = 1, heights = c(3, 1))

  list("Fig. 4g–i | Sequential splenic responses" = splenic_response_composite)
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
  cat("Fig. 4g-i | Splenic responses after sequential administration\nDescriptive statistics. Logarithmic panels summarize log10 measurements for display; their plotted centres differ from arithmetic means.\n")
  show("summary", "Arithmetic group means and variability", c("readout", "group", "n", "mean", "sd", "sem"))
  show("display_summary", "Values used for the plotted centres and error bars", c("readout", "group", "axis_scale", "n_display", "centre", "lower", "upper"))
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
    paths <- c(frames[basename(frames) == "Fig_04G-I_sequential_splenic_responses.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_04G-I_sequential_splenic_responses.R"], join(getwd(), "Fig_04G-I_sequential_splenic_responses.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_04G-I_sequential_splenic_responses.R"))], 1L)
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
      run <- join(base, "Fig_04G-I_sequential_splenic_responses", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
      needed <- c("dplyr", "ggbeeswarm", "ggplot2", "patchwork", "scales")
      missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
      if (length(missing)) stop("Install these plotting packages: ", paste(missing, collapse = ", "))
      if (interactive() && identical(Sys.getenv("RSTUDIO"), "1") && !identical(names(grDevices::dev.cur()), "RStudioGD"))
        grDevices::dev.new(noRStudioGD = FALSE)
      measure <- !interactive() && grDevices::dev.cur() == 1L
      if (measure) grDevices::pdf(NULL)
      plots <- tryCatch(suppressPackageStartupMessages(plot_results(analysis_result)), finally = if (measure) grDevices::dev.off())
      for (i in seq_along(plots)) {
        size <- analysis_result$sizes[[names(plots)[i]]]
        if (is.null(size)) size <- c(12, 6)
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
