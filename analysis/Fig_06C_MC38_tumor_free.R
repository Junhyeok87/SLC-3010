# Fig. 6c | Tumor-free MC38 mice at each assessment
# Count mice with a recorded tumor volume of zero on days 6, 12, 21 and 29, then divide
# by the number with a recorded volume on that day within each genotype.
# This panel pools the individual mice; it does not average experiment-specific percentages.
# Missing measurements are reported and never classified as tumor-free.
#
# Input: put Source_Data.xlsx beside this script or in the repository's data/ folder.
# Numerical packages: digest, readxl.
# Plotting packages: ggplot2.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$percentages
#   analysis_result$tables$individual_assessments
#   analysis_result$output_dir   # Folder containing Fig_06C_MC38_tumor_free_percentages.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Calculate the endpoint from mouse-level volumes, independently of the percentage summary.
analyze_data <- function(workbook, output_dir) {
  schedule <- data.frame(Cycle = 1:4, EvaluationDay = c(6, 12, 21, 29))
  strains <- c("C57BL/6", "Triple hIL-2")
  groups <- paste(strains, "SLC-3010", sep = "_")
  raw <- as.data.frame(readxl::read_excel(workbook, sheet = "Figure 6", range = "A7:K101",
                                         col_types = "text", .name_repair = "minimal"),
                       stringsAsFactors = FALSE)
  required <- c("Subject", "Group", paste("Day", schedule$EvaluationDay))
  names(raw) <- trimws(names(raw))
  if (anyDuplicated(names(raw)) || !all(required %in% names(raw)))
    stop("Expected distinct Subject, Group and assessment-day columns in the MC38 volume table.", call. = FALSE)
  blank <- function(x) is.na(x) | !nzchar(trimws(as.character(x)))
  raw <- raw[!Reduce(`&`, lapply(raw, blank)), , drop = FALSE]
  raw$Subject <- trimws(raw$Subject)
  raw$Group <- trimws(raw$Group)
  if (!nrow(raw) || any(blank(raw$Subject)) || any(blank(raw$Group)))
    stop("Every occupied MC38 row must identify the mouse and treatment group.", call. = FALSE)
  allowed <- c(groups, paste(strains, "PBS", sep = "_"))
  if (any(!raw$Group %in% allowed))
    stop("Unrecognized MC38 treatment group: ", paste(unique(raw$Group[!raw$Group %in% allowed]), collapse = ", "), call. = FALSE)
  if (anyDuplicated(raw[c("Group", "Subject")]))
    stop("Duplicate mouse identifier within an MC38 treatment group.", call. = FALSE)

  # The reported endpoint concerns the SLC-3010-treated cohorts. PBS cohorts are
  # not included in either the numerator or denominator of these proportions.
  raw <- raw[raw$Group %in% groups, , drop = FALSE]
  if (!all(groups %in% raw$Group)) stop("Both SLC-3010-treated MC38 genotypes are required.", call. = FALSE)
  raw$Strain <- strains[match(raw$Group, groups)]
  assessments <- do.call(rbind, lapply(seq_len(nrow(schedule)), function(i) {
    column <- paste("Day", schedule$EvaluationDay[i])
    text <- trimws(raw[[column]])
    missing <- is.na(text) | text %in% c("", "NA", "NaN", "N/A")
    volume <- suppressWarnings(as.numeric(text))
    invalid <- !missing & (!is.finite(volume) | volume < 0)
    if (any(invalid)) stop("Invalid MC38 volume on day ", schedule$EvaluationDay[i],
                           ", mouse ", paste(raw$Subject[invalid], collapse = ", "), call. = FALSE)
    volume[missing] <- NA_real_
    data.frame(Tumor = "MC38", Subject = raw$Subject, Group = raw$Group, Strain = raw$Strain,
               Cycle = schedule$Cycle[i], EvaluationDay = schedule$EvaluationDay[i],
               TumorVolume_mm3 = volume, Evaluated = !missing,
               TumorFree = ifelse(missing, NA, volume == 0), stringsAsFactors = FALSE)
  }))
  assessments <- assessments[order(assessments$Cycle, match(assessments$Strain, strains),
                                    assessments$Subject), , drop = FALSE]
  rownames(assessments) <- NULL
  percentages <- do.call(rbind, lapply(seq_len(nrow(schedule)), function(i) {
    do.call(rbind, lapply(strains, function(strain) {
      x <- assessments[assessments$Cycle == schedule$Cycle[i] & assessments$Strain == strain, ]
      n_evaluated <- sum(x$Evaluated)
      if (!n_evaluated) stop("No observed MC38 volumes for ", strain, " on day ",
                            schedule$EvaluationDay[i], ". A missing assessment is not 0% tumor-free.", call. = FALSE)
      n_free <- sum(x$TumorFree, na.rm = TRUE)
      pct <- 100 * n_free / n_evaluated
      data.frame(Tumor = "MC38", Strain = strain, Cycle = schedule$Cycle[i],
                 EvaluationDay = schedule$EvaluationDay[i], RegisteredMice = nrow(x),
                 EvaluatedMice = n_evaluated, MissingMice = sum(!x$Evaluated),
                 TumorFreeMice = n_free, TumorFree = pct, TumorBearing = 100 - pct,
                 stringsAsFactors = FALSE)
    }))
  }))
  rownames(percentages) <- NULL
  utils::write.csv(percentages, file.path(output_dir, "Fig_06C_MC38_tumor_free_percentages.csv"), row.names = FALSE)
  utils::write.csv(assessments, file.path(output_dir, "Fig_06C_MC38_individual_assessments.csv"), row.names = FALSE, na = "")
  utils::write.csv(schedule, file.path(output_dir, "Fig_06C_MC38_assessment_schedule.csv"), row.names = FALSE)
  list(tables = list(percentages = percentages, individual_assessments = assessments,
                     assessment_schedule = schedule))
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
  fig6c_bar_data <- function(data) {
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
  fig6c_make_plot <- function(d) {
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
          y = "SLC-3010 treatment cycle", title = "MC38") + g$guides(fill = g$guide_legend(ncol = 1, byrow = TRUE, 
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

  bars <- fig6c_bar_data(result$tables$percentages)
  plot <- fig6c_make_plot(bars)
  plot <- plot + ggplot2::labs(title = "Figure 6c | MC38 tumor-free mice",
    caption = paste("Tumor-free: recorded volume = 0 on days 6, 12, 21 and 29.",
                    "Percentages use pooled mouse counts within each genotype.", sep = "\n")) +
    ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0, size = 8))
  list("Figure 6c | MC38 tumor-free mice" = plot)
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
  cat("Fig. 6c | Tumor-free MC38 mice at each assessment\nBoth genotype cohorts received SLC-3010. Percentages use individual mice with observed measurements; missing measurements are not classified as tumor-free.\n")
  show("percentages", "Tumor-free counts, denominators and percentages", c("Strain", "Cycle", "EvaluationDay", "RegisteredMice", "EvaluatedMice", "MissingMice", "TumorFreeMice", "TumorFree"))
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
    paths <- c(frames[basename(frames) == "Fig_06C_MC38_tumor_free.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_06C_MC38_tumor_free.R"], join(getwd(), "Fig_06C_MC38_tumor_free.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_06C_MC38_tumor_free.R"))], 1L)
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
      run <- join(base, "Fig_06C_MC38_tumor_free", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
