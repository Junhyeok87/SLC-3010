# Fig. 6f | EO771 tumor-free percentages across experiments
# A mouse was tumor-free if its recorded tumor volume was zero on that cycle's evaluation day.
# Calculate the percentage separately in each experiment and average the three percentages
# with equal weight. The figure labels show that mean rounded to one decimal place.
#
# Input: put Source_Data.xlsx beside this script or in the repository's data/ folder.
# Numerical packages: digest, readxl.
# Plotting package: ggplot2.
# In RStudio, open this file and click Source to save the CSVs and then display the plot.
# To see the results, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$experiment_percentages  # Numerator, denominator and % in each experiment
#   analysis_result$tables$reported_percentages    # Exact mean and the percentage shown in Fig. 6f
#   analysis_result$output_dir   # Folder containing Fig_06F_EO771_tumor_free_percentages.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

analyze_data <- function(workbook, output_dir) {
  # Columns A:F are the experimental counts and their identifiers. Column G in Excel
  # displays the percentage for convenience; calculate it here directly from E and F.
  raw <- readxl::read_excel(workbook, sheet = "Figure 6", range = "A378:G402",
                            col_types = c(rep("text", 6), "skip"), .name_repair = "minimal")
  required <- c("Experiment", "Strain", "Cycle", "Evaluation day", "Tumor-free mice", "Evaluated mice")
  if (!identical(names(raw), required)) {
    stop("Fig. 6f requires the experiment-level mouse counts in Figure 6!A378:G402. ",
         "Select the matching Source_Data.xlsx accompanying this code release. ",
         "Expected columns: ", paste(required, collapse = ", "), ".", call. = FALSE)
  }

  as_whole_number <- function(column) {
    values <- suppressWarnings(as.numeric(raw[[column]]))
    if (any(!is.finite(values)) || any(values != floor(values))) {
      stop("Expected whole numbers in the Fig. 6f column '", column, "'.", call. = FALSE)
    }
    values
  }
  trials <- data.frame(
    Experiment = as_whole_number("Experiment"),
    Strain = trimws(raw$Strain),
    Cycle = as_whole_number("Cycle"),
    EvaluationDay = as_whole_number("Evaluation day"),
    TumorFreeMice = as_whole_number("Tumor-free mice"),
    EvaluatedMice = as_whole_number("Evaluated mice"),
    stringsAsFactors = FALSE
  )
  strains <- c("C57BL/6", "Triple hIL-2")
  evaluation_days <- c(6, 12, 21, 29)
  if (any(!trials$Experiment %in% 1:3) || any(!trials$Strain %in% strains) ||
      any(!trials$Cycle %in% seq_along(evaluation_days))) {
    stop("Fig. 6f requires experiments 1-3 and cycles 1-4 for C57BL/6 and Triple hIL-2.", call. = FALSE)
  }
  if (anyDuplicated(trials[c("Experiment", "Strain", "Cycle")]) ||
      nrow(trials) != 3L * length(strains) * length(evaluation_days)) {
    stop("Expected exactly one row per experiment, strain and cycle in Fig. 6f.", call. = FALSE)
  }
  if (any(trials$EvaluationDay != evaluation_days[trials$Cycle])) {
    stop("The four EO771 cycles were assessed on days 6, 12, 21 and 29, respectively.", call. = FALSE)
  }
  if (any(trials$EvaluatedMice <= 0 | trials$TumorFreeMice < 0 |
          trials$TumorFreeMice > trials$EvaluatedMice)) {
    stop("Check Fig. 6f counts: the denominator must be positive and the numerator between zero and the denominator.", call. = FALSE)
  }

  trials$TumorFreePct <- 100 * trials$TumorFreeMice / trials$EvaluatedMice
  trials <- trials[order(trials$Cycle, match(trials$Strain, strains), trials$Experiment), ]
  rownames(trials) <- NULL

  # Each of the three independent experiments contributes one percentage at each cycle.
  means <- stats::aggregate(TumorFreePct ~ Strain + Cycle + EvaluationDay, data = trials, FUN = mean)
  summary <- data.frame(
    Tumor = "EO771", Strain = means$Strain, Cycle = means$Cycle,
    TumorFree = round(means$TumorFreePct, 1), TumorBearing = 100 - round(means$TumorFreePct, 1),
    EvaluationDay = means$EvaluationDay, ExperimentCount = length(unique(trials$Experiment)),
    MeanTumorFreePct = means$TumorFreePct, stringsAsFactors = FALSE
  )
  summary <- summary[order(summary$Cycle, match(summary$Strain, strains)), ]
  rownames(summary) <- NULL
  utils::write.csv(trials, file.path(output_dir, "Fig_06F_EO771_tumor_free_by_experiment.csv"), row.names = FALSE)
  utils::write.csv(summary, file.path(output_dir, "Fig_06F_EO771_tumor_free_percentages.csv"), row.names = FALSE)
  list(tables = list(reported_percentages = summary, experiment_percentages = trials))
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result) {
  BASE_FAMILY <- "sans"
  AXIS_TEXT_PT <- 9
  AXIS_TITLE_PT <- 10
  PERCENT_LABEL_PT <- 6
  LEGEND_TEXT_PT <- 8.5
  TITLE_PT <- 7
  BAR_HEIGHT <- 0.31
  PAIR_OFFSET <- 0.19
  LABEL_INSIDE_MIN_PCT <- 25
  LABEL_PADDING_PCT <- 2
  FILL_COLORS <- c(`WT (C57BL/6), tumor-free` = "#4C78A8", `Triple hIL-2, tumor-free` = "#B91C1C", `Tumor-bearing` = "#D2D2D2")
  fig6f_bar_data <- function(data) {
      d <- data
      d$bar_y <- d$Cycle + ifelse(d$Strain == "C57BL/6", PAIR_OFFSET, -PAIR_OFFSET)
      d$ymin <- d$bar_y - BAR_HEIGHT/2
      d$ymax <- d$bar_y + BAR_HEIGHT/2
      d$free_key <- ifelse(d$Strain == "C57BL/6", "WT (C57BL/6), tumor-free", "Triple hIL-2, tumor-free")
      d$label <- ifelse(abs(d$TumorFree - round(d$TumorFree)) < 0.05, paste0(round(d$TumorFree), "%"), 
          paste0(formatC(d$TumorFree, format = "f", digits = 1), "%"))
      inside <- d$TumorFree >= LABEL_INSIDE_MIN_PCT
      d$label_x <- ifelse(inside, d$TumorFree - LABEL_PADDING_PCT, d$TumorFree + LABEL_PADDING_PCT)
      d$label_hjust <- ifelse(inside, 1, 0)
      d$label_color <- ifelse(inside, "#FFFFFF", "#222222")
      d
  }
  fig6f_make_plot <- function(d) {
      g <- asNamespace("ggplot2")
      p <- g$ggplot() + g$geom_hline(yintercept = c(1.5, 2.5, 3.5), color = "#E5E5E5", linewidth = 0.2) + 
          g$geom_rect(data = d, g$aes(xmin = TumorFree, xmax = 100, ymin = ymin, ymax = ymax, fill = "Tumor-bearing"), 
              color = NA, show.legend = TRUE) + g$geom_rect(data = d, g$aes(xmin = 0, xmax = TumorFree, 
          ymin = ymin, ymax = ymax, fill = free_key), color = NA, show.legend = TRUE) + g$geom_text(data = d, 
          g$aes(x = label_x, y = bar_y, label = label, hjust = label_hjust, color = label_color), family = BASE_FAMILY, 
          size = PERCENT_LABEL_PT/(72.27/25.4), vjust = 0.5, show.legend = FALSE) + g$scale_fill_manual(values = FILL_COLORS, 
          limits = names(FILL_COLORS), breaks = names(FILL_COLORS), drop = FALSE, name = NULL) + g$scale_color_identity() + 
          g$scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), expand = g$expansion(mult = 0)) + 
          g$scale_y_continuous(limits = c(0.5, 4.5), breaks = 1:4, expand = g$expansion(mult = 0)) + g$labs(x = "Mice (%)", 
          y = "SLC-3010 treatment cycle", title = "EO771") + g$guides(fill = g$guide_legend(ncol = 1, byrow = TRUE, 
          override.aes = list(color = NA))) + g$theme_bw(base_family = BASE_FAMILY, base_size = AXIS_TEXT_PT) + 
          g$theme(plot.background = g$element_rect(fill = "white", color = NA), panel.background = g$element_rect(fill = "white", 
              color = NA), panel.border = g$element_rect(fill = NA, color = "#444444", linewidth = 0.28), 
              panel.grid = g$element_blank(), axis.ticks = g$element_line(color = "#222222", linewidth = 0.26), 
              axis.ticks.length = grid::unit(1, "pt"), axis.text = g$element_text(size = AXIS_TEXT_PT, 
                  color = "#222222"), axis.text.x = g$element_text(margin = g$margin(t = 1.5)), axis.text.y = g$element_text(margin = g$margin(r = 1.5)), 
              axis.title = g$element_text(size = AXIS_TITLE_PT, color = "#111111"), axis.title.x = g$element_text(margin = g$margin(t = 2)), 
              axis.title.y = g$element_text(margin = g$margin(r = 2)), plot.title = g$element_text(size = TITLE_PT, 
                  face = "bold", hjust = 0.5, margin = g$margin(b = 2)), legend.position = "top", legend.justification = "left", 
              legend.direction = "vertical", legend.background = g$element_rect(fill = "white", color = NA), 
              legend.key = g$element_rect(fill = "white", color = NA), legend.key.width = grid::unit(6, 
                  "pt"), legend.key.height = grid::unit(5.8, "pt"), legend.text = g$element_text(size = LEGEND_TEXT_PT, 
                  color = "#222222", margin = g$margin(l = 1.5)), legend.spacing.y = grid::unit(0, "pt"), 
              legend.box.spacing = grid::unit(1.5, "pt"), legend.margin = g$margin(0, 0, 1.5, 0), plot.margin = g$margin(2, 
                  8, 2, 2))
      p
  }

  bars <- fig6f_bar_data(result$tables$reported_percentages)
  plot <- fig6f_make_plot(bars)
  plot <- plot + ggplot2::labs(title = "Figure 6f | EO771 tumor-free mice",
    caption = "Mean of the tumor-free percentages from three independent experiments, weighted equally.\nA mouse was tumor-free if its recorded tumor volume was zero on the evaluation day.") +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, size = 8))
  list("Figure 6f | EO771 tumor-free mice" = plot)
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
  cat("Fig. 6f | EO771 tumor-free percentages across experiments\nThe figure reports the equally weighted mean of the three experiment-specific percentages, not a pooled mouse proportion.\n")
  show("reported_percentages", "Equally weighted experiment means (%)", c("Strain", "Cycle", "EvaluationDay", "ExperimentCount", "MeanTumorFreePct"))
  show("experiment_percentages", "Experiment-specific numerators, denominators and percentages", c("Experiment", "Strain", "Cycle", "EvaluationDay", "TumorFreeMice", "EvaluatedMice", "TumorFreePct"))
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
    paths <- c(frames[basename(frames) == "Fig_06F_EO771_tumor_free.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_06F_EO771_tumor_free.R"], join(getwd(), "Fig_06F_EO771_tumor_free.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_06F_EO771_tumor_free.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply one workbook path or directory, or use Source in RStudio.", call. = FALSE)
    path <- path.expand(if (length(args)) args[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR"))
    book <- "Source_Data.xlsx"
    sheets <- c("Figure 6")
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
      run <- join(base, "Fig_06F_EO771_tumor_free", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
