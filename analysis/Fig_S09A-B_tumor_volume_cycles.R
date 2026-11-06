# Supplementary Fig. 9a-b | Tumor volumes by treatment cycle
# Within each tumor model and cycle, compare C57BL/6 and Triple hIL-2 mouse volumes with a
# two-sided, unpaired Welch t-test.
# The estimated difference is C57BL/6 minus Triple hIL-2. P values are nominal for each model
# and cycle.
#
# Input: put Source_Data.xlsx beside this script or in the repository's data/ folder.
# Numerical packages: digest, readxl.
# Plotting packages: ggplot2, scales.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$comparisons
#   analysis_result$output_dir   # Folder containing Fig_S09A-B_comparisons.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Compare the two treated genotypes separately at each assessment cycle.
# Each Welch contrast is C57BL/6 minus Triple hIL-2; P values are unadjusted.
# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  groups <- c("C57BL/6_SLC-3010", "Triple hIL-2_SLC-3010")
  ranges <- c(MC38 = "A7:F54", EO771 = "A60:F107")
  read_model <- function(model) {
    x <- as.data.frame(readxl::read_excel(workbook, sheet = "Supplementary Figure 9",
                                          range = ranges[[model]], col_types = "text", .name_repair = "minimal"))
    visits <- paste("Cycle", 1:4)
    if (!all(c("Subject", "Group", visits) %in% names(x)))
      stop("Missing Subject, Group, or Cycle columns for ", model, call. = FALSE)
    empty <- Reduce(`&`, lapply(x, function(z) is.na(z) | !nzchar(trimws(z))))
    x <- x[!empty, , drop = FALSE]
    x$Subject <- trimws(x$Subject)
    x$Group <- trimws(x$Group)
    if (!nrow(x) || anyNA(x$Subject) || any(!nzchar(x$Subject)) ||
        any(!x$Group %in% groups) || anyDuplicated(x[c("Group", "Subject")]))
      stop("Check mouse identifiers and treatment groups for ", model, call. = FALSE)
    # Preserve the original row order, cycling through the four observations per subject.
    cells <- as.vector(t(as.matrix(x[visits])))
    absent <- is.na(cells) | trimws(cells) %in% c("", "NA", "NaN")
    volume <- suppressWarnings(as.numeric(cells))
    if (any(!absent & (!is.finite(volume) | volume < 0)))
      stop("Invalid tumor volume for ", model, call. = FALSE)
    volume[absent] <- NA_real_
    data.frame(Model = model, Subject = rep(x$Subject, each = 4L), Group = rep(x$Group, each = 4L),
               Cycle = rep(1:4, nrow(x)), tumor_vol = volume, stringsAsFactors = FALSE)
  }
  assessments <- do.call(rbind, lapply(names(ranges), read_model))
  observed <- assessments[!is.na(assessments$tumor_vol), , drop = FALSE]

  compare_cycle <- function(model, cycle) {
    rows <- assessments[assessments$Model == model & assessments$Cycle == cycle, , drop = FALSE]
    x <- rows$tumor_vol[rows$Group == groups[1L] & !is.na(rows$tumor_vol)]
    y <- rows$tumor_vol[rows$Group == groups[2L] & !is.na(rows$tumor_vol)]
    fit <- if (min(length(x), length(y)) >= 2L) {
      tryCatch(stats::t.test(x, y, alternative = "two.sided", paired = FALSE, var.equal = FALSE),
               error = function(e) NULL)
    } else NULL
    mean_x <- if (length(x)) mean(x) else NA_real_
    mean_y <- if (length(y)) mean(y) else NA_real_
    p <- if (is.null(fit)) NA_real_ else fit$p.value
    status <- if (min(length(x), length(y)) < 2L) "Insufficient observations" else if (!is.finite(p))
      "Test not estimable" else if (p >= 0.05) "Not significant" else "Significant"
    data.frame(Model = model, Cycle = cycle, group1 = groups[1], group2 = groups[2],
               n_C57BL6 = length(x), n_Triple_hIL2 = length(y), mean_C57BL6 = mean_x, mean_Triple_hIL2 = mean_y,
               mean_difference = mean_x - mean_y,
               ci_low = if (is.null(fit)) NA_real_ else fit$conf.int[1L],
               ci_high = if (is.null(fit)) NA_real_ else fit$conf.int[2L],
               t = if (is.null(fit)) NA_real_ else unname(fit$statistic),
               df = if (is.null(fit)) NA_real_ else unname(fit$parameter), p_value = p, status = status,
               method = "Two-sided unpaired Welch t-test", adjustment = "none", stringsAsFactors = FALSE)
  }
  comparisons <- do.call(rbind, lapply(names(ranges), function(model) {
    do.call(rbind, lapply(1:4, function(cycle) compare_cycle(model, cycle)))
  }))
  coverage <- do.call(rbind, lapply(names(ranges), function(model) {
    do.call(rbind, lapply(groups, function(group) {
      do.call(rbind, lapply(1:4, function(cycle) {
        v <- assessments$tumor_vol[assessments$Model == model & assessments$Group == group & assessments$Cycle == cycle]
        data.frame(Model = model, Group = group, Cycle = cycle, n_registered = length(v),
                   n_observed = sum(!is.na(v)), n_missing = sum(is.na(v)))
      }))
    }))
  }))
  utils::write.csv(observed, file.path(output_dir, "Fig_S09A-B_individual_data.csv"), row.names = FALSE)
  utils::write.csv(comparisons, file.path(output_dir, "Fig_S09A-B_comparisons.csv"), row.names = FALSE)
  utils::write.csv(comparisons[c("Model", "Cycle", "n_C57BL6", "n_Triple_hIL2", "p_value", "status")],
                   file.path(output_dir, "Fig_S09A-B_test_status.csv"), row.names = FALSE)
  utils::write.csv(coverage, file.path(output_dir, "Fig_S09A-B_cycle_coverage.csv"), row.names = FALSE)
  list(tables = list(comparisons = comparisons, individual_measurements = observed, coverage = coverage))
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result, display = FALSE) {
  if (utils::packageVersion("ggplot2") < "4.0.0")
    stop("The cycle plot uses ggplot2 >= 4.0.0.", call. = FALSE)
  group_labels <- c("C57BL/6_SLC-3010" = "C57BL/6 + SLC-3010",
                    "Triple hIL-2_SLC-3010" = "Triple hIL-2 + SLC-3010")
  pal <- c("C57BL/6 + SLC-3010" = "#1F4CFF", "Triple hIL-2 + SLC-3010" = "#E31A1C")
  plot_df <- result$tables$individual_measurements
  plot_df$Model <- factor(plot_df$Model, levels = c("MC38", "EO771"))
  plot_df$Group_label <- factor(unname(group_labels[plot_df$Group]), levels = unname(group_labels))
  plot_df$x_num <- plot_df$Cycle
  plot_df$Cycle <- factor(plot_df$Cycle, levels = 1:4)
  stat_df <- result$tables$comparisons
  stat_df$x_num <- stat_df$Cycle
  pval <- stat_df$p_value
  stat_df$plot_label <- ifelse(!is.finite(pval), "NE", ifelse(pval < 1e-4, "****",
    ifelse(pval < 0.001, "***", ifelse(pval < 0.005, "**", ifelse(pval < 0.05, "*", "ns")))))
  stat_df$y_label <- vapply(seq_len(nrow(stat_df)), function(i) {
    m <- plot_df[as.character(plot_df$Model) == stat_df$Model[i], ]
    observed <- m$tumor_vol[m$x_num == stat_df$x_num[i]]
    if (!length(observed)) return(NA_real_)
    max(observed) + 0.1 * max(m$tumor_vol)
  }, numeric(1))
  stat_df$Model <- factor(stat_df$Model, levels = c("MC38", "EO771"))
  dodge_width <- 0.62
  jitter_width <- 0.11
  box_width <- 0.46
  x_min <- min(plot_df$x_num, na.rm = TRUE)
  x_max <- max(plot_df$x_num, na.rm = TRUE)
  x_left_limit <- x_min - 0.4
  x_right_limit <- x_max + dodge_width/4 + box_width/2
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x_num, y = tumor_vol, fill = Group_label, color = Group_label)) + ggplot2::geom_boxplot(ggplot2::aes(group = interaction(Cycle, 
      Group_label)), width = box_width, linewidth = 0.42, alpha = 0.18, outlier.shape = NA, position = ggplot2::position_dodge(width = dodge_width)) + 
      ggplot2::geom_point(ggplot2::aes(group = Group_label), position = ggplot2::position_jitterdodge(dodge.width = dodge_width, jitter.width = jitter_width, 
          jitter.height = 0, seed = 3010), size = 1.55, stroke = 0.22, alpha = 0.88) + ggplot2::stat_summary(ggplot2::aes(group = interaction(Cycle, 
      Group_label)), fun = median, geom = "crossbar", width = 0.34, middle_gp = list(linewidth = 0), linewidth = 0.45, color = "#2B2B2B", 
      position = ggplot2::position_dodge(width = dodge_width), show.legend = FALSE) + ggplot2::geom_text(data = stat_df, 
      ggplot2::aes(x = x_num, y = y_label, label = plot_label), inherit.aes = FALSE, size = if (display) 11 / ggplot2::.pt else 1.9, fontface = "bold", 
      color = "#111111", vjust = 0) + ggplot2::facet_wrap(~Model, nrow = 1, scales = "free_y") + ggplot2::scale_x_continuous(breaks = sort(unique(plot_df$x_num)), 
      labels = sort(unique(as.character(plot_df$Cycle))), limits = c(x_left_limit, x_right_limit), expand = ggplot2::expansion(mult = c(0, 
          0), add = c(0, 0))) + ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 1, big.mark = ","), 
      expand = ggplot2::expansion(mult = c(0.02, 0.18))) + ggplot2::scale_fill_manual(values = pal, name = NULL) + ggplot2::scale_color_manual(values = pal, 
      name = NULL) + ggplot2::labs(x = "SLC-3010 treatment cycle", y = expression("Tumor volume (mm"^3 * ")")) + 
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 0.35, size = 2), nrow = 1, byrow = TRUE), 
          color = "none") + ggplot2::coord_cartesian(clip = "off") + ggplot2::theme_classic(base_family = "sans", base_size = 6.6) + 
      ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA), panel.background = ggplot2::element_rect(fill = "white", 
          color = NA), panel.border = ggplot2::element_rect(fill = NA, color = "#2B2B2B", linewidth = 0.42), panel.grid.major.y = ggplot2::element_line(color = "#E6E6E6", 
          linewidth = 0.25), panel.grid.minor = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "#F0F0F0", 
          color = "#2B2B2B", linewidth = 0.36), strip.text = ggplot2::element_text(size = 6.5, face = "bold", color = "#222222", 
          margin = ggplot2::margin(1, 0, 1, 0)), axis.title.x = ggplot2::element_text(size = 7, face = "bold", margin = ggplot2::margin(t = 3.5)), 
          axis.title.y = ggplot2::element_text(size = 7, face = "bold", margin = ggplot2::margin(r = 3.5)), axis.text = ggplot2::element_text(size = 6.2, 
              color = "#222222"), axis.line = ggplot2::element_line(color = "#2B2B2B", linewidth = 0.38), axis.ticks = ggplot2::element_line(color = "#2B2B2B", 
              linewidth = 0.34), legend.position = "top", legend.justification = "right", legend.direction = "horizontal", 
          legend.box.just = "right", legend.background = ggplot2::element_blank(), legend.box.background = ggplot2::element_blank(), 
          legend.box.margin = ggplot2::margin(0, 0, -2, 0), legend.margin = ggplot2::margin(0, 0, 0, 0), legend.key.height = grid::unit(3.2, 
              "mm"), legend.key.width = grid::unit(4.8, "mm"), legend.text = ggplot2::element_text(size = 6.2, color = "#111111"), 
          panel.spacing.x = grid::unit(4, "mm"), plot.margin = ggplot2::margin(t = 1.2, r = 3.2, b = 2.2, l = 2.2, unit = "mm"))

  p <- p + ggplot2::labs(title = "Supplementary Figure 9a/b | Tumor volumes by treatment cycle",
    caption = "Two-sided unpaired Welch tests. ns: P ≥ 0.05. NE: comparison could not be estimated.\n* P < 0.05, ** P < 0.005, *** P < 0.001, **** P < 0.0001.") +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, size = 8))
  list("Supplementary Figure 9a/b | Tumor volumes by treatment cycle" = p)
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
  cat("Supplementary Fig. 9a-b | Tumor volumes by treatment cycle\nUnpaired, two-sided Welch t-tests. Difference: C57BL/6 minus Triple hIL-2. P values are nominal within each model and cycle.\n")
  show("comparisons", "Genotype differences and 95% confidence intervals", c("Model", "Cycle", "n_C57BL6", "n_Triple_hIL2", "mean_difference", "ci_low", "ci_high"))
  show("comparisons", "Welch t-tests by tumor model and cycle", c("Model", "Cycle", "t", "df", "p_value", "status"))
  show("coverage", "Registered, observed and missing mice", c("Model", "Group", "Cycle", "n_registered", "n_observed", "n_missing"))
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
    paths <- c(frames[basename(frames) == "Fig_S09A-B_tumor_volume_cycles.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_S09A-B_tumor_volume_cycles.R"], join(getwd(), "Fig_S09A-B_tumor_volume_cycles.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_S09A-B_tumor_volume_cycles.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply one workbook path or directory, or use Source in RStudio.", call. = FALSE)
    path <- path.expand(if (length(args)) args[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR"))
    book <- "Source_Data.xlsx"
    sheets <- c("Supplementary Figure 9")
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
      run <- join(base, "Fig_S09A-B_tumor_volume_cycles", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
      needed <- c("ggplot2", "scales")
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
        if (is.null(size)) size <- c(11, 6)
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
