# Fig. 5a-b | Vehicle-normalized blood Treg responses
# Summarize the supplied vehicle-normalized measurements by genotype, treatment and day using
# the mean and SD.
# The workbook supplies normalized values. Absolute Treg frequencies cannot be recovered from
# this table alone.
#
# Input: put Source_Data.xlsx beside this script or in the repository's data/ folder.
# Numerical packages: digest, dplyr, readxl.
# Plotting packages: dplyr, ggplot2, patchwork.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$summary
#   analysis_result$output_dir   # Folder containing Fig_05A-B_summary.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Individual vehicle-normalized Treg values; arithmetic mean +/- SD, fixed treatment groups and dosing markers.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  suppressPackageStartupMessages(library(dplyr))

  treatment_levels <- c("PBS", "hnTCB2", "SLC-3010")
  strain_levels <- c("C57BL/6", "Triple hIL-2")
  readout_keep <- "Relative Treg percentage"

  raw_data <- readxl::read_excel(workbook, sheet = "Figure 5", range = "A7:G381", col_types = c("numeric", "text", "text", "numeric", "text", "text", "numeric"), .name_repair = "minimal")

  names(raw_data) <- trimws(names(raw_data))
  required_cols <- c("mouse_id", "strain", "group", "day", "treatment", "readout", "value")
  missing_cols <- setdiff(required_cols, names(raw_data))
  if (length(missing_cols) > 0) {
      stop("Input file is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  measurements <- raw_data %>% dplyr::select(dplyr::all_of(required_cols)) %>% dplyr::mutate(mouse_id = as.character(mouse_id), 
      strain = as.character(strain), group = as.character(group), day = as.numeric(day), treatment = as.character(treatment), 
      readout = as.character(readout), value = as.numeric(value)) %>% dplyr::filter(!is.na(strain), !is.na(day), 
      !is.na(treatment), !is.na(readout), !is.na(value)) %>% dplyr::filter(readout == readout_keep) %>% 
      dplyr::filter(strain %in% strain_levels, treatment %in% treatment_levels) %>% dplyr::mutate(strain = factor(strain, 
      levels = strain_levels), treatment = factor(treatment, levels = treatment_levels))
  if (nrow(measurements) == 0) {
      stop("No rows matched the expected readout, strain and treatment. In Source_Data.xlsx, open Data_index and check the rows listed for Fig_05A-B_blood_Treg.")
  }
  missing_strains <- setdiff(strain_levels, unique(as.character(measurements$strain)))
  if (length(missing_strains) > 0) {
      warning("These strain(s) were not found: ", paste(missing_strains, collapse = ", "))
  }
  missing_treatments <- setdiff(treatment_levels, unique(as.character(measurements$treatment)))
  if (length(missing_treatments) > 0) {
      warning("These treatment(s) were not found: ", paste(missing_treatments, collapse = ", "))
  }

  summary_data <- measurements %>% dplyr::group_by(strain, treatment, day) %>% dplyr::summarise(n = dplyr::n(), 
      mean = mean(value, na.rm = TRUE), sd = stats::sd(value, na.rm = TRUE), sem = sd/sqrt(n), error = sd, 
      ymin = pmax(mean - error, 0), ymax = mean + error, .groups = "drop") %>% dplyr::mutate(strain = factor(as.character(strain), 
      levels = strain_levels), treatment = factor(as.character(treatment), levels = treatment_levels))

  utils::write.csv(measurements, file.path(output_dir, "Fig_05A-B_individual_data.csv"), row.names = FALSE)
  utils::write.csv(summary_data, file.path(output_dir, "Fig_05A-B_summary.csv"), row.names = FALSE)

  list(tables = list(summary = summary_data, observations = measurements))
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result) {
  suppressPackageStartupMessages(library(dplyr))

  plot_data <- result$tables$observations
  summary_data <- result$tables$summary
  treatment_levels <- levels(plot_data$treatment)
  strain_levels <- levels(plot_data$strain)

  base_family <- "Arial"
  admin_days <- c(0, 7, 14, 21)
  admin_marker_label <- "SLC-3010/hnTCB2 administration"
  admin_marker_color <- "#FF1717"
  admin_marker_shape <- 17
  admin_marker_size <- 1.65
  admin_marker_y_fraction <- -0.08

  plot_labels <- c(`C57BL/6` = "C57BL/6", `Triple hIL-2` = "Triple hIL-2")
  treatment_colors <- c(PBS = "#8A8A8A", hnTCB2 = "#3346FF", `SLC-3010` = "#FF1717")
  treatment_shapes <- c(PBS = 16, hnTCB2 = 17, `SLC-3010` = 15)
  axis_label_x <- expression("Days after" ~ 1^"st" ~ "SLC-3010 or hnTCB2 administration")
  axis_label_y <- "Relative Treg percentage\n(Normalized to vehicle ctrl.)"
  line_width <- 0.42
  point_size <- 1.55
  errorbar_width <- 0.32
  errorbar_line <- 0.38

  all_days <- sort(unique(plot_data$day))
  x_major_breaks <- c(0, 7, 14, 21, 28, 35)
  x_axis_min <- min(c(all_days, admin_days, x_major_breaks), na.rm = TRUE)
  x_axis_max <- max(c(all_days, admin_days, x_major_breaks), na.rm = TRUE)

  nice_upper <- function(x) {
      x <- max(x, na.rm = TRUE)
      if (!is.finite(x) || x <= 0) 
          return(1)
      raw <- x * 1.05
      pow <- 10^floor(log10(raw))
      candidates <- c(1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10) * pow
      out <- min(candidates[candidates >= raw])
      if (!is.finite(out)) 
          out <- ceiling(raw)
      out
  }
  y_limit_tbl <- summary_data %>% dplyr::group_by(strain) %>% dplyr::summarise(max_summary = max(ymax, 
      na.rm = TRUE), .groups = "drop") %>% dplyr::left_join(plot_data %>% dplyr::group_by(strain) %>% dplyr::summarise(max_raw = max(value, 
      na.rm = TRUE), .groups = "drop"), by = "strain") %>% dplyr::mutate(y_max = vapply(pmax(max_summary, 
      max_raw), nice_upper, numeric(1)))
  get_y_max <- function(strain_name) {
      y <- y_limit_tbl$y_max[as.character(y_limit_tbl$strain) == strain_name]
      if (length(y) == 0 || !is.finite(y)) 
          stop("Could not determine y-axis maximum for ", strain_name)
      y[[1]]
  }
  key_levels <- c(treatment_levels, admin_marker_label)
  key_colors <- c(treatment_colors, setNames(admin_marker_color, admin_marker_label))
  key_shapes <- c(treatment_shapes, setNames(admin_marker_shape, admin_marker_label))
  make_panel <- function(strain_name) {
    values <- summary_data[as.character(summary_data$strain) == strain_name, ]
    upper <- get_y_max(strain_name)
    doses <- data.frame(day = admin_days, y = admin_marker_y_fraction * upper, key = admin_marker_label)
    ggplot2::ggplot(values, ggplot2::aes(day, mean, colour = treatment, shape = treatment)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax),
                             width = errorbar_width, linewidth = errorbar_line, show.legend = FALSE) +
      ggplot2::geom_line(ggplot2::aes(group = treatment), linewidth = line_width, show.legend = FALSE) +
      ggplot2::geom_point(size = point_size, show.legend = TRUE) +
      ggplot2::geom_point(data = doses, ggplot2::aes(day, y, colour = key, shape = key),
                          inherit.aes = FALSE, size = admin_marker_size, show.legend = FALSE) +
      ggplot2::scale_colour_manual(name = NULL, values = key_colors, limits = key_levels, drop = FALSE) +
      ggplot2::scale_shape_manual(name = NULL, values = key_shapes, limits = key_levels, drop = FALSE) +
      ggplot2::scale_x_continuous(breaks = x_major_breaks, limits = c(x_axis_min - 0.65, x_axis_max + 0.55)) +
      ggplot2::coord_cartesian(ylim = c(admin_marker_y_fraction * upper, upper), clip = "off") +
      ggplot2::labs(title = strain_name, x = axis_label_x, y = axis_label_y) +
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(hjust = 0.5))
  }
  panels <- patchwork::wrap_plots(lapply(strain_levels, make_panel), nrow = 1, guides = "collect") +
    patchwork::plot_annotation(caption = "Mean ± SD. Treg frequencies are normalized to genotype-matched vehicle controls.")
  panels <- panels & ggplot2::theme(legend.position = "bottom")

  list("Fig. 5a–b | Blood Treg responses" = panels)
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
  cat("Fig. 5a-b | Vehicle-normalized blood Treg responses\nDescriptive statistics of the supplied vehicle-normalized measurements; the figure uses mean and SD.\n")
  show("summary", "Normalized Treg measurements by genotype, treatment and day", c("strain", "treatment", "day", "n", "mean", "sd", "sem"))
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
    paths <- c(frames[basename(frames) == "Fig_05A-B_blood_Treg.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_05A-B_blood_Treg.R"], join(getwd(), "Fig_05A-B_blood_Treg.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_05A-B_blood_Treg.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply one workbook path or directory, or use Source in RStudio.", call. = FALSE)
    path <- path.expand(if (length(args)) args[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR"))
    book <- "Source_Data.xlsx"
    sheets <- c("Figure 5")
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
      run <- join(base, "Fig_05A-B_blood_Treg", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
      needed <- c("dplyr", "ggplot2", "patchwork")
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
