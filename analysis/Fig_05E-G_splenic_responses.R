# Fig. 5e-g | Genotype comparisons of splenic responses
# For each outcome, compare the two genotypes on days 4, 11, 18 and 25 with separate
# two-sided, unpaired Welch t-tests.
# The estimated difference is Triple hIL-2 minus C57BL/6. P values are nominal for each day.
#
# Input: put Source_Data.xlsx beside this script or in the repository's data/ folder.
# Numerical packages: digest, dplyr, readxl.
# Plotting packages: dplyr, ggplot2, patchwork.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$comparisons
#   analysis_result$output_dir   # Folder containing Fig_05E-G_test_status.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Splenic MP CD8 T cells, Tregs and ratios; arithmetic mean +/- SD and unpaired Welch tests on days 4, 11, 18 and 25.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  suppressPackageStartupMessages(library(dplyr))

  stat_days <- c(4, 11, 18, 25)
  stat_compare_strains <- c("C57BL/6", "Triple hIL-2")

  target_strains <- c("C57BL/6", "Triple hIL-2")
  readout_levels <- c("MP_CD8T", "Treg", "MP_CD8T_Treg_ratio")

  raw_data <- dplyr::bind_rows(
      readxl::read_excel(workbook, sheet = "Figure 5", range = "A387:F507", col_types = c("numeric", "text", "numeric", "text", "text", "numeric"), .name_repair = "minimal"),
      readxl::read_excel(workbook, sheet = "Figure 5", range = "A513:F633", col_types = c("numeric", "text", "numeric", "text", "text", "numeric"), .name_repair = "minimal"),
      readxl::read_excel(workbook, sheet = "Figure 5", range = "A639:F759", col_types = c("numeric", "text", "numeric", "text", "text", "numeric"), .name_repair = "minimal")
  )

  names(raw_data) <- trimws(names(raw_data))
  required_cols <- c("mouse_id", "strain", "day", "treatment", "readout", "value")
  missing_cols <- setdiff(required_cols, names(raw_data))
  if (length(missing_cols) > 0) {
      stop("Input file is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  normalize_strain <- function(x) {
      x0 <- trimws(as.character(x))
      x_lower <- tolower(x0)
      out <- ifelse(grepl("c57|b6|wt|wild", x_lower), "C57BL/6", x0)
      out <- ifelse(grepl("triple|hil[- ]?2|hil2|human.*il[- ]?2", x_lower), "Triple hIL-2", out)
      out
  }
  normalize_readout <- function(x) {
      x0 <- trimws(as.character(x))
      x_lower <- tolower(gsub("[[:space:]/-]+", "_", x0))
      out <- ifelse(grepl("mp.*cd8", x_lower) & grepl("treg|ratio", x_lower), "MP_CD8T_Treg_ratio", x0)
      out <- ifelse(grepl("treg", x_lower) & !grepl("ratio", x_lower), "Treg", out)
      out <- ifelse(grepl("mp.*cd8", x_lower) & !grepl("treg|ratio", x_lower), "MP_CD8T", out)
      out
  }
  measurements <- raw_data %>% dplyr::select(dplyr::all_of(required_cols)) %>% dplyr::mutate(mouse_id = as.character(mouse_id), 
      strain = normalize_strain(strain), day = as.numeric(day), treatment = trimws(as.character(treatment)), 
      readout = normalize_readout(readout), value = as.numeric(value)) %>% dplyr::filter(strain %in% target_strains, 
      readout %in% readout_levels, !is.na(day), !is.na(value)) %>% dplyr::mutate(strain = factor(strain, 
      levels = target_strains), readout = factor(readout, levels = readout_levels))

  selected <- !is.na(measurements$treatment) & grepl("^SLC[-_ ]?3010$", measurements$treatment, ignore.case = TRUE)
  if (!any(selected)) stop("This analysis requires observations of SLC-3010 treatment in both strains.")
  measurements <- measurements[selected, , drop = FALSE]
  if (anyDuplicated(measurements[c("mouse_id", "strain", "day", "readout")]))
    stop("Duplicate mouse, strain, day and readout keys in SLC-3010 observations.")

  if (nrow(measurements) == 0) stop("No usable rows remain after cleaning/filtering data.")
  if (!all(target_strains %in% unique(as.character(measurements$strain)))) {
      stop("Both C57BL/6 and Triple hIL-2 data are required; single-strain outputs are intentionally not generated.")
  }
  summary_data <- measurements %>% dplyr::group_by(strain, readout, day) %>% dplyr::summarise(n = dplyr::n(), 
      mean = mean(value, na.rm = TRUE), sd = stats::sd(value, na.rm = TRUE), sem = sd/sqrt(n), error = sd, 
      ymin = pmax(mean - error, 0), ymax = mean + error, .groups = "drop") %>% dplyr::mutate(sd = ifelse(is.na(sd), 
      0, sd), sem = ifelse(is.na(sem), 0, sem), error = ifelse(is.na(error), 0, error), ymin = ifelse(is.na(ymin), 
      mean, ymin), ymax = ifelse(is.na(ymax), mean, ymax))
  p_to_label <- function(p) {
      dplyr::case_when(is.na(p) ~ NA_character_, p < 1e-04 ~ "****", p < 0.001 ~ "***", p < 0.005 ~ "**", 
          p < 0.05 ~ "*", TRUE ~ "ns")
  }

  compare_strains <- function(d) {
    wt <- d$value[d$strain == stat_compare_strains[1L]]
    triple <- d$value[d$strain == stat_compare_strains[2L]]
    test <- tryCatch(stats::t.test(triple, wt, paired = FALSE, var.equal = FALSE,
                                  alternative = "two.sided", conf.level = 0.95), error = identity)
    ok <- !inherits(test, "error")
    data.frame(n_C57BL6 = length(wt), n_Triple_hIL2 = length(triple),
               mean_C57BL6 = mean(wt), mean_Triple_hIL2 = mean(triple),
               sd_C57BL6 = sd(wt), sd_Triple_hIL2 = sd(triple),
               effect_direction = "Triple hIL-2 minus C57BL/6",
               mean_difference = mean(triple) - mean(wt), confidence_level = 0.95,
               ci_low = if (ok) test$conf.int[1L] else NA_real_,
               ci_high = if (ok) test$conf.int[2L] else NA_real_,
               t_statistic = if (ok) unname(test$statistic) else NA_real_,
               df = if (ok) unname(test$parameter) else NA_real_,
               p_value = if (ok) test$p.value else NA_real_,
               p_adjustment = "none", status_detail = if (ok) "Estimated" else conditionMessage(test))
  }
  stat_results <- measurements %>%
    filter(day %in% stat_days) %>%
    group_by(readout, day) %>% group_modify(~ compare_strains(.x)) %>% ungroup() %>%
    mutate(p_label = p_to_label(p_value), test = "unpaired Welch t-test",
           comparison = "Triple hIL-2 vs C57BL/6")

  utils::write.csv(measurements, file.path(output_dir, "Fig_05E-G_individual_data.csv"), row.names = FALSE)
  utils::write.csv(summary_data, file.path(output_dir, "Fig_05E-G_summary.csv"), row.names = FALSE)
  utils::write.csv(stat_results, file.path(output_dir, "Fig_05E-G_comparisons.csv"), row.names = FALSE)
  by_day <- transform(stat_results, test_status = ifelse(is.na(p_value), "Not estimable",
                       ifelse(p_label == "ns", "Not significant", "Significant")))
  utils::write.csv(by_day, file.path(output_dir, "Fig_05E-G_test_status.csv"), row.names = FALSE)

  list(tables = list(comparisons = by_day, summary = summary_data, observations = measurements))
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result) {
  suppressPackageStartupMessages(library(dplyr))

  stat_results <- result$tables$comparisons
  summary_data <- result$tables$summary %>% mutate(legend_item = as.character(strain))
  readout_levels <- levels(summary_data$readout)
  stat_days <- sort(unique(stat_results$day))

  stat_text_size <- 2.1
  stat_text_vjust <- 0.5
  stat_y_pad_fraction <- 0.075
  stat_y_cap_fraction <- 0.955
  admin_days <- c(0, 7, 14, 21)
  admin_color <- "#FF1717"
  admin_marker_y_fraction <- -0.055

  legend_levels <- c("C57BL/6", "Triple hIL-2", "SLC-3010 administration")
  readout_ylabels <- c(MP_CD8T = "MP CD8 T cells (x10⁶ cells)", Treg = "Treg (x10⁶ cells)", MP_CD8T_Treg_ratio = "MP CD8 T/Treg ratio")
  y_max_by_readout <- c(MP_CD8T = 300, Treg = 20, MP_CD8T_Treg_ratio = 80)
  y_breaks <- list(MP_CD8T = c(0, 100, 200, 300), Treg = c(0, 5, 10, 15, 20), MP_CD8T_Treg_ratio = c(0, 
      20, 40, 60, 80))
  legend_colors <- c(`C57BL/6` = "#525252", `Triple hIL-2` = "#D7301F", `SLC-3010 administration` = admin_color)
  legend_shapes <- c(`C57BL/6` = 16, `Triple hIL-2` = 15, `SLC-3010 administration` = 17)
  legend_linetypes <- c(`C57BL/6` = "solid", `Triple hIL-2` = "solid", `SLC-3010 administration` = "blank")

  stat_annotation_data <- summary_data %>% dplyr::filter(day %in% stat_days) %>% dplyr::group_by(readout, 
      day) %>% dplyr::summarise(y_raw = max(ymax, na.rm = TRUE), .groups = "drop") %>% dplyr::left_join(stat_results, 
      by = c("readout", "day")) %>% dplyr::mutate(y_max = y_max_by_readout[as.character(readout)], y = pmin(y_raw + 
      stat_y_pad_fraction * y_max, stat_y_cap_fraction * y_max), p_label = ifelse(is.na(p_label), "NE", p_label)) %>% 
      dplyr::filter(nzchar(p_label))
  make_panel <- function(readout_name) {
    df <- summary_data[as.character(summary_data$readout) == readout_name, ]
    upper <- unname(y_max_by_readout[readout_name])
    doses <- data.frame(day = admin_days, y = admin_marker_y_fraction * upper,
                        key = "SLC-3010 administration")
    annotations <- stat_annotation_data[as.character(stat_annotation_data$readout) == readout_name, ]
    ggplot2::ggplot(df, ggplot2::aes(day, mean, colour = legend_item, shape = legend_item)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax), width = 0.32,
                             linewidth = 0.31, show.legend = FALSE) +
      ggplot2::geom_line(ggplot2::aes(group = legend_item), linewidth = 0.5, show.legend = FALSE) +
      ggplot2::geom_point(size = 1.7, show.legend = TRUE) +
      ggplot2::geom_point(data = doses, ggplot2::aes(day, y, colour = key, shape = key),
                          inherit.aes = FALSE, size = 2, show.legend = FALSE) +
      ggplot2::geom_text(data = annotations, ggplot2::aes(day, y, label = p_label),
                         inherit.aes = FALSE, size = 3, fontface = "bold") +
      ggplot2::scale_color_manual(name = NULL, values = legend_colors, limits = legend_levels, drop = FALSE) +
      ggplot2::scale_shape_manual(name = NULL, values = legend_shapes, limits = legend_levels, drop = FALSE) +
      ggplot2::scale_x_continuous(breaks = c(0, 7, 14, 21, 28, 35), limits = c(-0.6, 35.6)) +
      ggplot2::scale_y_continuous(breaks = y_breaks[[readout_name]]) +
      ggplot2::coord_cartesian(ylim = c(admin_marker_y_fraction * upper, upper), clip = "off") +
      ggplot2::labs(x = expression("Days after" ~ 1^"st" ~ "SLC-3010 administration"),
                    y = readout_ylabels[[readout_name]]) +
      ggplot2::theme_classic(base_size = 10)
  }
  panels <- patchwork::wrap_plots(lapply(readout_levels, make_panel), nrow = 1, guides = "collect") +
    patchwork::plot_annotation(caption = "Mean ± SD. Between-strain Welch tests on days 4, 11, 18 and 25. ns: not significant; NE: not estimable.")
  panels <- panels & ggplot2::theme(legend.position = "bottom")

  list("Fig. 5e–g | Splenic responses" = panels)
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
  cat("Fig. 5e-g | Genotype comparisons of splenic responses\nUnpaired, two-sided Welch t-tests. Difference: Triple hIL-2 minus C57BL/6. P values are nominal for each readout and day.\n")
  show("comparisons", "Genotype differences and 95% confidence intervals", c("readout", "day", "n_C57BL6", "n_Triple_hIL2", "mean_difference", "ci_low", "ci_high"))
  show("comparisons", "Welch t-tests by readout and day", c("readout", "day", "t_statistic", "df", "p_value", "status_detail"))
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
    paths <- c(frames[basename(frames) == "Fig_05E-G_splenic_responses.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_05E-G_splenic_responses.R"], join(getwd(), "Fig_05E-G_splenic_responses.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_05E-G_splenic_responses.R"))], 1L)
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
      run <- join(base, "Fig_05E-G_splenic_responses", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
        if (is.null(size)) size <- c(11, 5)
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
