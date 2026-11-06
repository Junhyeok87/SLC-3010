# Fig. 3a | MO5 growth and day-21 transfer comparisons
# Summarize observed tumor volumes by group and day, then compare mice with and without
# transfer within Beta IL-2 and SLC-3010 using two-sided, unpaired equal-variance t tests.
# The day-21 comparisons report unadjusted P values. The growth summary lists the number of
# observed mice at each time.
#
# Input: put Source_Data.xlsx beside this script or in the repository's data/ folder.
# Numerical packages: digest, readxl.
# Plotting packages: ggplot2.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$comparisons
#   analysis_result$output_dir   # Folder containing Fig_03A_day21_ttest_stats.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Figure 3a: MO5 growth, with or without OT-1 transfer.
# Reads the individual-animal Source Data table.
# Arithmetic means and SEM are computed from observed volumes at each day;
# blank source cells remain missing and are not imputed.
# Day-21 brackets use the unadjusted P value from the two-sided unpaired Student
# t-test comparing mice with and without OT-1 transfer within each treatment.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  raw <- as.data.frame(readxl::read_excel(workbook, sheet = "Figure 3", range = "A8:G776", .name_repair = "minimal"),
                       stringsAsFactors = FALSE)
  required <- c("Group", "OT1_transfer", "Treatment", "Replicate", "Mouse_ID",
                "Day", "Tumor_volume_mm3")
  if (!identical(names(raw), required)) stop("Unexpected Fig. 3a source columns.")
  if (nrow(raw) != 768L) stop("Expected all 96 within-group replicate positions at all eight days.")
  group_levels <- paste0("G", 1:6)
  days <- seq(0, 21, 3)
  group_titles <- c("w/o OT-1+Vehicle", "w/o OT-1+Beta IL-2", "w/o OT-1+SLC-3010",
                    "w/ OT-1+Vehicle", "w/ OT-1+Beta IL-2", "w OT-1+SLC-3010")

  if (!all(raw$Group %in% group_levels) || !all(raw$Day %in% days))
    stop("Unrecognized group or day in the source table.")
  group_index <- match(raw$Group, group_levels)
  if (any(raw$Treatment != rep(c("Vehicle", "Beta IL-2", "SLC-3010"), 2)[group_index]) ||
      any(raw$OT1_transfer != rep(c("No", "Yes"), each = 3)[group_index]))
    stop("Source group identities do not match the treatment and OT-1 transfer definitions.")
  key <- paste(raw$Group, raw$Replicate, raw$Day, sep = "|")
  if (anyDuplicated(key)) stop("Duplicate animal/day source observations.")
  if (any(!is.na(raw$Tumor_volume_mm3) &
          (!is.finite(raw$Tumor_volume_mm3) | raw$Tumor_volume_mm3 < 0)))
    stop("Tumor volumes must be nonnegative finite numbers or missing.")
  summary <- do.call(rbind, lapply(seq_along(group_levels), function(g) {
    do.call(rbind, lapply(days, function(day) {
      rows <- raw[raw$Group == group_levels[g] & raw$Day == day, , drop = FALSE]
      if (nrow(rows) != 16L || !identical(sort(as.integer(rows$Replicate)), 1:16))
        stop("Expected every within-group replicate position, including missing measurements, at each group/day.")
      values <- rows$Tumor_volume_mm3[!is.na(rows$Tumor_volume_mm3)]
      n <- length(values)
      data.frame(Group = group_levels[g], Group_label = group_titles[g], Day = day,
        n = n, Mean = if (n) mean(values) else NA_real_,
        SD = if (n > 1L) stats::sd(values) else NA_real_,
        SEM = if (n > 1L) stats::sd(values) / sqrt(n) else NA_real_)
    }))
  }))
  write.csv(summary, file.path(output_dir,
    "Fig_03A_MO5_tumor_growth_summary.csv"), row.names = FALSE, na = "")

  comparisons <- data.frame(Treatment = c("Beta IL-2", "SLC-3010"),
    Group1 = c("G2", "G3"), Group2 = c("G5", "G6"), Day = 21,
    stringsAsFactors = FALSE)
  test_stats <- do.call(rbind, lapply(seq_len(nrow(comparisons)), function(i) {
    a <- raw$Tumor_volume_mm3[raw$Group == comparisons$Group1[i] & raw$Day == 21]
    b <- raw$Tumor_volume_mm3[raw$Group == comparisons$Group2[i] & raw$Day == 21]
    a <- a[!is.na(a)]; b <- b[!is.na(b)]
    test <- stats::t.test(a, b, alternative = "two.sided", paired = FALSE,
                          var.equal = TRUE)
    p <- as.numeric(test$p.value)
    label <- if (p < 0.0001) "****" else if (p < 0.001) "***" else
      if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
    data.frame(comparisons[i, ], OT1_Group1 = "No", OT1_Group2 = "Yes",
      n1 = length(a), n2 = length(b), Mean1 = mean(a), Mean2 = mean(b),
      effect_direction = "without minus with OT-1", mean_difference = mean(a) - mean(b),
      confidence_level = 0.95, ci_low = test$conf.int[1L], ci_high = test$conf.int[2L],
      t = as.numeric(test$statistic), df = as.numeric(test$parameter),
      p_value = p, label = label, test_method = "Student's t-test",
      alternative = "two.sided", paired = FALSE, var_equal = TRUE,
      p_adjustment = "none", stringsAsFactors = FALSE)
  }))
  write.csv(test_stats, file.path(output_dir,
    "Fig_03A_day21_ttest_stats.csv"), row.names = FALSE)

  list(tables = list(comparisons = test_stats, summary = summary, observations = raw))
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result) {
  summary <- result$tables$summary
  test_stats <- result$tables$comparisons
  group_levels <- unique(summary$Group)
  days <- sort(unique(summary$Day))

  group_colors <- c("#000000", "#0000C0", "#FF0000", "#000000", "#0000C0", "#A00000")

  # Use the same means, SEM and day-21 comparisons in a standard R plot.
  key_labels <- paste(rep(c("Without OT-1", "With OT-1"), each = 3),
                      rep(c("Vehicle", "Beta IL-2", "SLC-3010"), 2), sep = " | ")
  display <- transform(summary, Group = factor(Group, levels = group_levels))
  brackets <- transform(test_stats, x = c(21.55, 22.45),
                        ymin = pmin(Mean1, Mean2), ymax = pmax(Mean1, Mean2))
  p <- ggplot2::ggplot(display, ggplot2::aes(Day, Mean, colour = Group, shape = Group)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Mean - SEM, ymax = Mean + SEM),
                           width = 0.35, linewidth = 0.4, na.rm = TRUE, show.legend = FALSE) +
    ggplot2::geom_line(linewidth = 0.55, na.rm = TRUE) +
    ggplot2::geom_point(size = 2, na.rm = TRUE) +
    ggplot2::scale_colour_manual(name = NULL, values = setNames(group_colors, group_levels),
                                 labels = key_labels) +
    ggplot2::scale_shape_manual(name = NULL, values = setNames(rep(c(1, 17), each = 3), group_levels),
                                labels = key_labels) +
    ggplot2::geom_segment(data = brackets, ggplot2::aes(x = x, xend = x, y = ymin, yend = ymax),
                           inherit.aes = FALSE, linewidth = 0.4) +
    ggplot2::geom_segment(data = brackets, ggplot2::aes(x = x - 0.22, xend = x, y = ymin, yend = ymin),
                           inherit.aes = FALSE, linewidth = 0.4) +
    ggplot2::geom_segment(data = brackets, ggplot2::aes(x = x - 0.22, xend = x, y = ymax, yend = ymax),
                           inherit.aes = FALSE, linewidth = 0.4) +
    ggplot2::geom_text(data = brackets, ggplot2::aes(x = x + 0.14, y = (ymin + ymax) / 2, label = label),
                        inherit.aes = FALSE, hjust = 0, size = 3) +
    ggplot2::scale_x_continuous(breaks = days, limits = c(0, 23.5)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 4000, 1000)) +
    ggplot2::coord_cartesian(ylim = c(0, max(4000, display$Mean + display$SEM, na.rm = TRUE))) +
    ggplot2::labs(x = "Days after OT-1 cell transfer", y = expression(Tumor~volume~(mm^3)),
                  caption = "Mean ± SEM. Day 21: without vs with OT-1 within each indicated treatment.") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 8),
                   plot.caption = ggplot2::element_text(hjust = 0)) +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 3, byrow = TRUE),
                    shape = ggplot2::guide_legend(ncol = 3, byrow = TRUE))

  list("Fig. 3a | MO5 tumor growth" = p)
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
  cat("Fig. 3a | MO5 growth and day-21 transfer comparisons\nUnpaired, two-sided Student t-tests with equal variances; P values are unadjusted. Difference: without minus with OT-1 transfer.\n")
  show("comparisons", "Day-21 mean differences and 95% confidence intervals", c("Treatment", "OT1_Group1", "OT1_Group2", "n1", "n2", "Mean1", "Mean2", "mean_difference", "ci_low", "ci_high"))
  show("comparisons", "Day-21 Student t-tests", c("Treatment", "t", "df", "p_value"))
  cat("\nComplete result tables: ", paste(names(result$tables), collapse = ", "), "\n", sep = "")
  cat("Results directory: ", result$output_dir, "\n", sep = "")
  invisible(result)
}

# Source calculates the results, prints the analysis summary and displays the plots.
if (isTRUE(getOption("SLC3010.autorun", TRUE))) {
  analysis_result <- figure_result <- NULL
  analysis_result <- local({
    needed <- c("digest", "readxl")
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
    paths <- c(frames[basename(frames) == "Fig_03A_MO5_tumor_growth.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_03A_MO5_tumor_growth.R"], join(getwd(), "Fig_03A_MO5_tumor_growth.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_03A_MO5_tumor_growth.R"))], 1L)
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
      run <- join(base, "Fig_03A_MO5_tumor_growth", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
      needed <- c("ggplot2")
      missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
      if (length(missing)) stop("Install these plotting packages: ", paste(missing, collapse = ", "))
      if (interactive() && identical(Sys.getenv("RSTUDIO"), "1") && !identical(names(grDevices::dev.cur()), "RStudioGD"))
        grDevices::dev.new(noRStudioGD = FALSE)
      measure <- !interactive() && grDevices::dev.cur() == 1L
      if (measure) grDevices::pdf(NULL)
      plots <- tryCatch(suppressPackageStartupMessages(plot_results(analysis_result)), finally = if (measure) grDevices::dev.off())
      for (i in seq_along(plots)) {
        size <- analysis_result$sizes[[names(plots)[i]]]
        if (is.null(size)) size <- c(10, 6)
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
