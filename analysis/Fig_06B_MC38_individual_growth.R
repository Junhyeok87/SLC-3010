# Fig. 6b | Individual MC38 tumor trajectories
# Organize mouse-level tumor volumes by genotype, treatment and study day and plot each
# mouse's observed measurements.
#
# Input: put Source_Data.xlsx beside this script or in the repository's data/ folder.
# Numerical packages: digest, readxl.
# Plotting packages: ggplot2, scales.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$individual_measurements
#   analysis_result$output_dir   # Folder containing Fig_06B_MC38_individual_growth_individual_data.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Observational unit: a mouse within its assigned group, measured on a study day.
# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  group_order <- c("Triple hIL-2 + PBS", "Triple hIL-2 + SLC-3010",
                   "C57BL/6 + PBS", "C57BL/6 + SLC-3010")
  raw <- as.data.frame(readxl::read_excel(workbook, sheet = "Figure 6", range = "A7:K101",
                                         col_types = "text", .name_repair = "minimal"))
  clean_text <- function(x) {
      x <- enc2utf8(as.character(x))
      x <- gsub("[   ]", " ", x, perl = TRUE)
      x <- gsub("[​‌‍﻿]", "", x, perl = TRUE)
      x <- gsub("[‐‑‒–—−]", "-", x, perl = TRUE)
      trimws(gsub("[[:space:]]+", " ", x))
  }
  names(raw) <- clean_text(names(raw))
  if (anyDuplicated(names(raw)) || !all(c("Subject", "Group") %in% names(raw)))
    stop("Expected distinct column names including Subject and Group.", call. = FALSE)
  visits <- setdiff(names(raw), c("Subject", "Group"))
  days <- suppressWarnings(as.numeric(sub("^Day ", "", visits)))
  if (!length(days) || any(!is.finite(days)) || anyDuplicated(days))
    stop("Expected a unique numeric study day for each measurement column.", call. = FALSE)
  ord <- order(days)
  visits <- visits[ord]
  days <- days[ord]
  blank <- function(x) is.na(x) | !nzchar(trimws(as.character(x)))
  raw <- raw[!Reduce(`&`, lapply(raw, blank)), , drop = FALSE]
  subject <- clean_text(raw$Subject)
  group_input <- clean_text(raw$Group)
  if (!nrow(raw) || any(blank(subject)) || any(blank(group_input)))
    stop("Every occupied data row must identify the mouse and treatment group.", call. = FALSE)
  group_key <- function(x) tolower(gsub("[[:space:]_+]", "", clean_text(x)))
  group <- group_order[match(group_key(group_input), group_key(group_order))]
  if (anyNA(group))
    stop("Unrecognized treatment group: ", paste(unique(group_input[is.na(group)]), collapse = ", "))
  if (anyDuplicated(data.frame(group, subject)))
    stop("Duplicate subject within a group; use the original mouse and experiment identifiers.", call. = FALSE)

  # Retain true zero volumes. Missing measurements contribute to the coverage
  # table but do not enter arithmetic summaries or individual trajectories.
  values <- lapply(seq_along(visits), function(i) {
    cell <- trimws(raw[[visits[i]]])
    absent <- is.na(cell) | cell %in% c("", "NA", "NaN")
    volume <- suppressWarnings(as.numeric(cell))
    invalid <- !absent & (!is.finite(volume) | volume < 0)
    if (any(invalid)) stop("Invalid tumor volume on day ", days[i], ", mouse ",
                           paste(subject[invalid], collapse = ", "), call. = FALSE)
    volume[absent] <- NA_real_
    data.frame(subject = subject, group_label = group, day = days[i], tumor_volume = volume,
               stringsAsFactors = FALSE)
  })
  all_visits <- do.call(rbind, values)
  measurements <- all_visits[!is.na(all_visits$tumor_volume), , drop = FALSE]
  if (length(setdiff(group_order, measurements$group_label)))
    stop("Each study group must contain at least one observed tumor volume.", call. = FALSE)
  measurements$strain <- ifelse(startsWith(measurements$group_label, "C57BL/6"), "C57BL/6", "Triple hIL-2")
  measurements$treatment <- ifelse(endsWith(measurements$group_label, "SLC-3010"), "SLC-3010", "PBS")
  measurements <- measurements[order(match(measurements$group_label, group_order),
                                      measurements$subject, measurements$day), ]
  rownames(measurements) <- NULL

  summary <- do.call(rbind, lapply(group_order, function(g) {
    do.call(rbind, lapply(days, function(day) {
      visit <- all_visits[all_visits$group_label == g & all_visits$day == day, ]
      x <- visit$tumor_volume[!is.na(visit$tumor_volume)]
      n <- length(x)
      s <- if (n >= 2L) stats::sd(x) else NA_real_
      data.frame(group_label = g, day = day, n_registered = nrow(visit), n = n,
                 n_missing = nrow(visit) - n, mean = if (n) mean(x) else NA_real_,
                 sd = s, sem = if (n >= 2L) s / sqrt(n) else NA_real_,
                 median = if (n) stats::median(x) else NA_real_)
    }))
  }))
  utils::write.csv(measurements, file.path(output_dir, "Fig_06B_MC38_individual_growth_individual_data.csv"), row.names = FALSE)
  utils::write.csv(summary, file.path(output_dir, "Fig_06B_MC38_individual_growth_summary.csv"), row.names = FALSE)
  list(tables = list(group_summary = summary, individual_measurements = measurements))
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result) {
  GROUP_LEVELS <- c("Triple hIL-2 + PBS", "Triple hIL-2 + SLC-3010", "C57BL/6 + PBS", "C57BL/6 + SLC-3010")
  FACET_LEVELS <- c("C57BL/6 + PBS", "C57BL/6 + SLC-3010", "Triple hIL-2 + PBS", "Triple hIL-2 + SLC-3010")
  GROUP_COLORS <- c(`Triple hIL-2 + PBS` = "#1F1F1F", `Triple hIL-2 + SLC-3010` = "#FF1717", `C57BL/6 + PBS` = "#1F1F1F", 
      `C57BL/6 + SLC-3010` = "#2457FF")
  GROUP_SHAPES <- c(`Triple hIL-2 + PBS` = 18, `Triple hIL-2 + SLC-3010` = 18, `C57BL/6 + PBS` = 16, `C57BL/6 + SLC-3010` = 16)
  GROUP_LINETYPES <- c(`Triple hIL-2 + PBS` = "22", `Triple hIL-2 + SLC-3010` = "solid", `C57BL/6 + PBS` = "solid", 
      `C57BL/6 + SLC-3010` = "solid")
  FACET_LABELS <- c(`C57BL/6 + PBS` = "C57BL/6\nPBS", `C57BL/6 + SLC-3010` = "C57BL/6\nSLC-3010", `Triple hIL-2 + PBS` = "Triple hIL-2\nPBS", 
      `Triple hIL-2 + SLC-3010` = "Triple hIL-2\nSLC-3010")
  data <- result$tables$individual_measurements
  summary <- result$tables$group_summary
  summary <- summary[summary$n > 0, , drop = FALSE]
  data$group_label <- factor(data$group_label, levels = GROUP_LEVELS)
  summary$group_label <- factor(summary$group_label, levels = GROUP_LEVELS)
  data$facet_label <- factor(data$group_label, levels = FACET_LEVELS)
  summary$facet_label <- factor(summary$group_label, levels = FACET_LEVELS)
  values_to_show <- data$tumor_volume
  upper <- max(4000, ceiling(max(values_to_show, na.rm = TRUE) / 1000) * 1000)

  lower <- 0
  admin_y <- upper * 0.025
  breaks <- pretty(c(lower, upper), n = 5)
  y_scale <- ggplot2::scale_y_continuous(breaks = breaks, labels = scales::label_comma())
  admin_days <- c(0, 7, 14, 21)
  administration <- data.frame(day = admin_days, y = admin_y, key = "SLC-3010 administration")
  x_breaks <- sort(unique(c(0, 7, 14, 21, 28)))
  x_limits <- range(c(data$day, x_breaks)) + c(-0.9, 0.9)

  # Repeat the study timeline across facets as in the article; the label
  # distinguishes scheduled dosing of treated mice from the PBS controls.
  admin_key <- "SLC-3010 dosing schedule (treated groups)"
  administration$key <- admin_key
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = day, y = tumor_volume, group = interaction(group_label, subject))) +
    ggplot2::geom_line(ggplot2::aes(colour = group_label, linetype = group_label),
                       linewidth = 0.32, alpha = 0.72, show.legend = FALSE) +
    ggplot2::geom_point(ggplot2::aes(colour = group_label, shape = group_label),
                        size = 0.5, alpha = 0.45, show.legend = FALSE) +
    ggplot2::geom_line(data = summary, ggplot2::aes(x = day, y = median, group = facet_label,
                       colour = group_label, linetype = group_label), inherit.aes = FALSE,
                       linewidth = 0.72, alpha = 0.96, show.legend = FALSE) +
    ggplot2::geom_point(data = administration, ggplot2::aes(x = day, y = y, shape = key),
                        inherit.aes = FALSE, colour = "#FF1717", size = 2.1, show.legend = TRUE) +
    ggplot2::facet_wrap(~facet_label, nrow = 1, drop = FALSE, labeller = ggplot2::as_labeller(FACET_LABELS)) +
    ggplot2::scale_colour_manual(NULL, values = GROUP_COLORS, guide = "none") +
    ggplot2::scale_linetype_manual(NULL, values = GROUP_LINETYPES, guide = "none") +
    ggplot2::scale_shape_manual(NULL, values = c(GROUP_SHAPES, stats::setNames(17, admin_key)),
                                limits = c(GROUP_LEVELS, admin_key), breaks = admin_key)

  plot <- plot + y_scale +
    ggplot2::coord_cartesian(xlim = x_limits, ylim = c(lower, upper)) +
    ggplot2::scale_x_continuous(breaks = x_breaks, minor_breaks = NULL) +
    ggplot2::labs(title = "Figure 6b | MC38", subtitle = "Individual mice (thin lines) and group medians (thick lines).",
                  x = "Days after 1st SLC-3010 administration",
                  y = expression("Tumor volume (mm"^3 * ")"),
                  caption = "Full observed range. The article used a display cap of 3,980 mm³.") +
    ggplot2::theme_classic(base_size = 10, base_family = "sans") +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 8.5),
                   plot.caption = ggplot2::element_text(hjust = 0, size = 8),
                   strip.background = ggplot2::element_rect(fill = "#F4F4F4", color = "#727272"),
                   strip.text = ggplot2::element_text(size = 9))

  list("Figure 6b | MC38" = plot)
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
  cat("Fig. 6b | Individual MC38 tumor trajectories\nIndividual trajectories remain in the complete measurement table; missing observations are counted explicitly.\n")
  show("group_summary", "Observed cohort size and tumor volume by day (cubic millimetres)", c("group_label", "day", "n_registered", "n", "n_missing", "mean", "sd", "median"))
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
    paths <- c(frames[basename(frames) == "Fig_06B_MC38_individual_growth.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_06B_MC38_individual_growth.R"], join(getwd(), "Fig_06B_MC38_individual_growth.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_06B_MC38_individual_growth.R"))], 1L)
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
      run <- join(base, "Fig_06B_MC38_individual_growth", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
