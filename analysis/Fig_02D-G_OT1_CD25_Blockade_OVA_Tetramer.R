# Fig. 2d-g | CD25 blockade and OVA-tetramer CD8 responses
# Summarize counts by treatment for each panel, using the observed mice in d/e and the total
# and tetramer-positive counts recorded together in each f/g source row.
# The corresponding e, f and g summaries are returned separately. These panels report
# descriptive statistics.
#
# Input: put Source_Data.xlsx beside this script or in the repository's data/ folder.
# Numerical packages: digest, readxl.
# Plotting packages: ggplot2, patchwork.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$panel_d_summary
#   analysis_result$output_dir   # Folder containing Fig_02D_summary.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Figure 2d-g: CD25-blockade experiment and OVA-tetramer CD8 T cell responses.
# Panels d/e retain median, type-7 quartiles, 1.5-IQR whiskers and every raw observation.
# Panels f/g retain arithmetic mean +/- SEM and every workbook row.
# No hypothesis tests or significance labels are added. Input values are already in 10^6 cells.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  table_ranges <- c(
    "fig02d_ot1_blockade" = "'Figure 2'!A49:C73",
    "fig02e_host_blockade" = "'Figure 2'!A79:C103",
    "fig02fg_ova_tetramer" = "'Figure 2'!A109:D126"
  )

  read_panel_rows <- function(x) {
      x <- as.data.frame(x, stringsAsFactors = FALSE)
      blank <- function(z) is.na(z) | !nzchar(trimws(z))
      empty <- apply(x, 1, function(row) all(blank(row)))
      x <- x[!empty, , drop = FALSE]
      if (!nrow(x) || any(vapply(x, function(z) any(blank(z)), logical(1))))
          stop("Missing required Source Data values; no partial rows are discarded.")
      x[] <- lapply(x, trimws)
      x
  }

  numeric_values <- function(z, label) {
      result <- suppressWarnings(as.numeric(z))
      if (any(!is.finite(result) | result < 0)) 
          stop("Invalid value(s) in ", label)
      result
  }

  summarise_groups <- function(data) {
      result <- lapply(split(data, data$treatment_order), function(x) data.frame(treatment_order = x$treatment_order[1], n = nrow(x), 
          mean = mean(x$value), sd = sd(x$value), sem = sd(x$value)/sqrt(nrow(x)), median = median(x$value), q1 = unname(quantile(x$value, 
              0.25)), q3 = unname(quantile(x$value, 0.75)), min = min(x$value), max = max(x$value)))
      do.call(rbind, result)
  }

  fg_source <- readxl::read_excel(workbook, sheet = "Figure 2", range = "A109:D126", col_types = rep("text", 4), .name_repair = "minimal")
  fg_raw <- read_panel_rows(fg_source[c("treatment", "peptide", "total_cd8_count_1e6",
      "ova_tetramer_cd8_count_1e6")])

  fg_treatments <- c("Vehicle", "Beta IL-2", "SLC-3010")

  fg_peptides <- c("Vehicle", "SIINFEKL")

  if (any(!fg_raw$treatment %in% fg_treatments) || any(!fg_raw$peptide %in% fg_peptides)) {
      stop("Unknown treatment or peptide in the Figure 2f/g Source Data block")
  }

  fg_raw$treatment_order <- match(fg_raw$treatment, fg_treatments) + 3 * (match(fg_raw$peptide, fg_peptides) - 1)

  fg_raw$family <- c("Vehicle", "BetaIL2", "SLC3010")[match(fg_raw$treatment, fg_treatments)]

  if (!setequal(fg_raw$treatment_order, 1:6)) stop("A f/g group is absent")

  fg_panels <- lapply(c(f = "total_cd8_count_1e6", g = "ova_tetramer_cd8_count_1e6"), function(column) {
      x <- fg_raw
      x$value <- numeric_values(x[[column]], column)
      x
  })

  fg_summary <- lapply(fg_panels, summarise_groups)

  de_ids <- c("fig02d_ot1_blockade", "fig02e_host_blockade")
  de_source <- lapply(de_ids, function(id) readxl::read_excel(workbook, range = table_ranges[[id]], col_types = rep("text", 3), .name_repair = "minimal"))
  de_rows <- Map(function(x, population) data.frame(mouse_id = x$mouse_id,
      treatment = x$treatment, readout = population, value = x$cell_count_1e6,
      stringsAsFactors = FALSE), de_source, c("Thy1.1+CD8", "Thy1.1-CD8"))
  de_raw <- read_panel_rows(do.call(rbind, de_rows))

  de_raw$treatment <- gsub("\\s+", "", de_raw$treatment)

  de_raw$treatment[de_raw$treatment == "PBS"] <- "Vehicle"

  de_raw$treatment[de_raw$treatment %in% c("OVA+SLC-3010+PC61", "OVA+SLC-3010+Agly-PC61")] <- "OVA+SLC-3010+AglyPC61"

  de_levels <- c("Vehicle", "OVA", "SLC-3010", "OVA+SLC-3010", "OVA+SLC-3010+AglyPC61")

  de_readouts <- c(d = "Thy1.1+CD8", e = "Thy1.1-CD8")

  if (any(!de_raw$treatment %in% de_levels) || any(!de_raw$readout %in% de_readouts)) {
      stop("Unexpected treatment or population in the Figure 2d/e Source Data blocks.")
  }

  if (anyDuplicated(paste(de_raw$mouse_id, de_raw$readout, sep = "|"))) {
      stop("Duplicate mouse/readout observations in d/e input")
  }

  de_raw$value <- numeric_values(de_raw$value, "d/e value")

  de_raw$treatment_order <- match(de_raw$treatment, de_levels)

  de_raw <- de_raw[order(match(de_raw$readout, de_readouts), de_raw$treatment_order, suppressWarnings(as.numeric(de_raw$mouse_id)), 
      de_raw$mouse_id), ]

  de_panels <- lapply(de_readouts, function(readout) {
      x <- de_raw[de_raw$readout == readout, , drop = FALSE]
      if (!setequal(x$treatment_order, 1:5)) 
          stop("A d/e treatment group is absent")
      x
  })

  if (!setequal(de_panels$d$mouse_id, de_panels$e$mouse_id)) stop("d/e mouse IDs differ")

  de_summary <- lapply(de_panels, summarise_groups)

  panel_order <- c("d", "e", "f", "g")
  all_panels <- c(de_panels, fg_panels)

  all_summaries <- c(de_summary, fg_summary)
  numeric_tables <- list()

  for (panel in panel_order) {
      source <- all_panels[[panel]]
      summary <- all_summaries[[panel]]
      summary$treatment <- vapply(summary$treatment_order, function(i) source$treatment[match(i, source$treatment_order)], 
          character(1))
      summary$panel <- panel
      if (panel %in% c("f", "g")) 
          summary$peptide <- vapply(summary$treatment_order, function(i) source$peptide[match(i, source$treatment_order)], 
              character(1))
      numeric_tables[[paste0("panel_", panel, "_summary")]] <- summary
      write.csv(summary, file.path(output_dir, paste0("Fig_02", toupper(panel), "_summary.csv")), row.names = FALSE, 
          fileEncoding = "UTF-8")
  }


  de_values <- de_raw[, c("mouse_id", "treatment", "readout", "value")]
  fg_values <- data.frame(source_row = seq_len(nrow(fg_raw)), treatment = fg_raw$treatment,
      peptide = fg_raw$peptide, total_cd8_count_1e6 = fg_panels$f$value,
      ova_tetramer_cd8_count_1e6 = fg_panels$g$value)
  # A source row links the two f/g measurements; it is not an animal identifier.
  write.csv(de_values, file.path(output_dir, "Fig_02D-E_individual_counts.csv"), row.names = FALSE)
  write.csv(fg_values, file.path(output_dir, "Fig_02F-G_individual_counts.csv"), row.names = FALSE)
  numeric_tables$blockade_measurements <- de_values
  numeric_tables$tetramer_measurements <- fg_values
  list(tables = numeric_tables, de_panels = de_panels, fg_panels = fg_panels,
       de_levels = de_levels, peptides = fg_peptides, fg_summary = fg_summary)
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result) {
  library(ggplot2)

  font_family <- "Arial"
  axis_title_pt <- 7
  tick_pt <- 6
  axis_line_mm <- 0.22
  point_size_mm <- 1.3
  point_stroke_mm <- 0.25

  family_colours <- c(Vehicle = "#4D4D4D", OVA = "#000000", BetaIL2 = "#3B5FD9", SLC3010 = "#E53935")

  fg_family_order <- c("Vehicle", "BetaIL2", "SLC3010")

  fg_labels <- c("Vehicle", "Beta IL-2", "SLC-3010")
  fg_shapes <- c(Vehicle = 21, BetaIL2 = 17, SLC3010 = 15)

  fg_fills <- c(Vehicle = "white", BetaIL2 = "#3B5FD9", SLC3010 = "#E53935")

  fg_x_positions <- c(1, 2, 3, 4.8, 5.8, 6.8)

  fg_x_limits <- c(0.4, 7.4)

  spread_offsets <- function(n, half_width = 0.13) if (n == 1L) 0 else seq(-half_width, half_width, length.out = n)


  fg_panels <- result$fg_panels
  fg_summary <- result$fg_summary
  fg_peptides <- result$peptides
  de_levels <- result$de_levels
  de_panels <- result$de_panels
  fg_panels <- lapply(fg_panels, function(x) {
      x$x_group <- fg_x_positions[x$treatment_order]
      x$x_point <- x$x_group + ave(x$treatment_order, x$treatment_order,
          FUN = function(g) spread_offsets(length(g), 0.23))
      x
  })
  de_panels <- lapply(de_panels, function(x) {
      x$x_point <- x$treatment_order + ave(x$treatment_order, x$treatment_order,
          FUN = function(g) spread_offsets(length(g), 0.18))
      x
  })

  de_labels <- paste0("G", seq_along(de_levels), " ",
      c("Vehicle", "OVA", "SLC-3010", "OVA + SLC-3010", "OVA + SLC-3010 + Agly PC61"))

  de_colours <- setNames(c("#4D4D4D", "#000000", "#E53935", "#E53935", "#7B3294"), de_levels)

  de_fills <- setNames(c("white", "#000000", "#E53935", "#E53935", "#7B3294"), de_levels)

  de_shapes <- setNames(c(21, 16, 15, 16, 18), de_levels)

  theme_figure <- function(vertical_x = FALSE) theme_classic(base_size = tick_pt, base_family = font_family) + theme(axis.title.x = element_blank(), 
      axis.title.y = element_text(size = axis_title_pt, margin = margin(r = 1.7)), axis.text = element_text(size = tick_pt, 
          colour = "black"), axis.text.y = element_text(margin = margin(r = 0.9)), axis.text.x = element_text(angle = if (vertical_x) 90 else 0, 
          hjust = if (vertical_x) 1 else 0.5, vjust = 0.5, margin = margin(t = 1.5)), axis.line = element_line(linewidth = axis_line_mm, 
          colour = "black"), axis.ticks = element_line(linewidth = axis_line_mm, colour = "black"), axis.ticks.length = unit(1.2, 
          "pt"), legend.position = "bottom", plot.title = element_blank(), plot.margin = margin(0, 0, 0, 0, unit = "mm"))

  make_fg <- function(panel, ymax, breaks, y_title) {
      x <- fg_panels[[panel]]
      s <- fg_summary[[panel]]
      s$x_group <- fg_x_positions[s$treatment_order]
      observed_max <- max(c(x$value, s$mean + s$sem), na.rm = TRUE)
      if (observed_max > ymax) {
          ymax <- observed_max
          breaks <- pretty(c(0, ymax))
      }
      ggplot(x) + geom_errorbar(data = s, aes(x = x_group, ymin = mean - sem, ymax = mean + sem), width = 0.24, linewidth = 0.25, 
          colour = "black") + geom_segment(data = s, aes(x = x_group - 0.23, xend = x_group + 0.23, y = mean, yend = mean), 
          linewidth = 0.32, colour = "black") + geom_point(aes(x = x_point, y = value, colour = family, fill = family, 
          shape = family), size = point_size_mm, stroke = point_stroke_mm) + scale_colour_manual(name = NULL, values = family_colours, breaks = fg_family_order, labels = fg_labels) +
          scale_fill_manual(name = NULL, values = fg_fills, breaks = fg_family_order, labels = fg_labels) +
          scale_shape_manual(name = NULL, values = fg_shapes, breaks = fg_family_order, labels = fg_labels) + scale_x_continuous(breaks = c(2, 5.8), labels = fg_peptides, limits = fg_x_limits, expand = expansion(mult = 0)) + 
          scale_y_continuous(breaks = breaks, limits = c(0, ymax), expand = expansion(mult = c(0, 0.025))) + labs(y = y_title) + 
          coord_cartesian(clip = "off") + theme_figure(FALSE)
  }

  fg_plots <- list(f = make_fg("f", 150, seq(0, 150, 50), "CD8⁺ cells (×10⁶)"), g = make_fg("g", 10, seq(0, 
      10, 2), "OVA-tet⁺CD8⁺ cells (×10⁶)"))

  de_ymax <- c(d = 25, e = 300)

  make_de <- function(panel, breaks, y_title) {
      x <- de_panels[[panel]]
      ymax <- max(de_ymax[[panel]], x$value)
      if (ymax > de_ymax[[panel]]) breaks <- pretty(c(0, ymax))
      ggplot(x) + geom_boxplot(aes(x = treatment_order, y = value, group = treatment, fill = treatment), orientation = "x", 
          width = 0.48, linewidth = 0.24, colour = "black", alpha = 0.16, outlier.shape = NA, coef = 1.5) + geom_point(aes(x = x_point, 
          y = value, colour = treatment, fill = treatment, shape = treatment), size = point_size_mm, stroke = point_stroke_mm) + 
          scale_colour_manual(name = NULL, values = de_colours, breaks = de_levels, labels = de_labels) +
          scale_fill_manual(name = NULL, values = de_fills, breaks = de_levels, labels = de_labels) +
          scale_shape_manual(name = NULL, values = de_shapes, breaks = de_levels, labels = de_labels) + scale_x_continuous(breaks = 1:5, labels = paste0("G", 
          1:5), limits = c(0.55, 5.55), expand = expansion(mult = 0)) + scale_y_continuous(breaks = breaks, limits = c(0, 
          ymax), expand = expansion(mult = c(0, 0.025))) + labs(y = y_title) + coord_cartesian(clip = "off") + 
          theme_figure(FALSE)
  }

  de_plots <- list(d = make_de("d", seq(0, 25, 5), "Thy1.1⁺CD8⁺ (×10⁶)"), e = make_de("e", seq(0, 300, 100), 
      "Thy1.1⁻CD8⁺ (×10⁶)"))

  de_figure <- patchwork::wrap_plots(de_plots, nrow = 1, guides = "collect") &
      theme(legend.position = "bottom")
  fg_figure <- patchwork::wrap_plots(fg_plots, nrow = 1, guides = "collect") &
      theme(legend.position = "bottom")

  list("Fig. 2d-e - CD25 blockade" = de_figure, "Fig. 2f-g - OVA tetramer" = fg_figure)
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
  cat("Fig. 2d-g | CD25 blockade and OVA-tetramer CD8 responses\nPanels d/e display medians and quartiles; f/g display means and SEM. These are descriptive analyses.\n")
  show("panel_d_summary", "Panel d: transferred CD8 counts (millions)", c("treatment", "n", "median", "q1", "q3", "mean", "sd"))
  show("panel_e_summary", "Panel e: host CD8 counts (millions)", c("treatment", "n", "median", "q1", "q3", "mean", "sd"))
  show("panel_f_summary", "Panel f: total CD8 counts (millions)", c("treatment", "peptide", "n", "mean", "sd", "sem"))
  show("panel_g_summary", "Panel g: OVA-tetramer-positive CD8 counts (millions)", c("treatment", "peptide", "n", "mean", "sd", "sem"))
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
    paths <- c(frames[basename(frames) == "Fig_02D-G_OT1_CD25_Blockade_OVA_Tetramer.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_02D-G_OT1_CD25_Blockade_OVA_Tetramer.R"], join(getwd(), "Fig_02D-G_OT1_CD25_Blockade_OVA_Tetramer.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_02D-G_OT1_CD25_Blockade_OVA_Tetramer.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply one workbook path or directory, or use Source in RStudio.", call. = FALSE)
    path <- path.expand(if (length(args)) args[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR"))
    book <- "Source_Data.xlsx"
    sheets <- c("Figure 2")
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
      run <- join(base, "Fig_02D-G_OT1_CD25_Blockade_OVA_Tetramer", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
      needed <- c("ggplot2", "patchwork")
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
