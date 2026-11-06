# Supplementary Fig. 5b | Tissue radiotracer uptake
# Summarize mouse measurements as mean percentage of injected dose per gram and SEM within
# each tracer, tissue and time point.
# The workbook supplies zero values at time zero. The first scan was at 15 minutes.
#
# Input: put Source_Data.xlsx beside this script or in the repository's data/ folder.
# Numerical packages: digest, dplyr, readr, readxl.
# Plotting packages: ggh4x, ggplot2, scales.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$uptake_summary
#   analysis_result$output_dir   # Folder containing Fig_S05B_mean_SEM.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Supplementary Figure 5b: serial in vivo PET/CT biodistribution by tissue, with individual values and mean +/- SEM.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  library(dplyr)
  organ_levels <- c("Heart", "Intestine", "Kidney", "Liver", "Lung", "Spleen", "Stomach", "UB")
  group_levels <- c("IL-2/[⁸⁹Zr]hnTCB2", "[⁸⁹Zr]hnTCB2")

  raw <- readxl::read_excel(workbook, sheet = "Supplementary Figure 5", range = "A8:G776", col_types = c("numeric", "text", "numeric", "text", "numeric", "numeric", "text"), .name_repair = "minimal")

  if (!"mode" %in% names(raw) || anyNA(raw$mode) || any(raw$mode != "in vivo")) {
    stop("Supplementary Figure 5b requires in vivo PET/CT measurements for every row.", call. = FALSE)
  }

  dat <- raw %>%
    transmute(
      group   = as.character(group),
      subject = factor(subject),
      organ   = factor(organ, levels = organ_levels),
      t       = as.numeric(t),
      conc    = as.numeric(conc)
    ) %>%
    mutate(
      group_label = dplyr::recode(
        group,
        "SLC-3010" = "IL-2/[⁸⁹Zr]hnTCB2",
        "hnTCB2"     = "[⁸⁹Zr]hnTCB2",
        .default   = group
      ),
      group_label = factor(group_label, levels = group_levels)
    ) %>%
    filter(!is.na(t), !is.na(conc))
  if (anyNA(dat[c("group_label", "organ", "subject")]))
    stop("Missing or unrecognized tracer, tissue or subject in uptake measurements.")
  if (anyDuplicated(dat[c("group", "subject", "organ", "t")]))
    stop("A mouse has duplicate uptake measurements for the same tissue and time.")

  sum_dat <- dat %>%
    group_by(group_label, organ, t) %>%
    summarise(
      mean_conc = mean(conc, na.rm = TRUE),
      sd_conc   = sd(conc, na.rm = TRUE),
      n         = dplyr::n(),
      se_conc   = sd_conc / sqrt(n),
      .groups   = "drop"
    ) %>%
    arrange(organ, group_label, t)

  readr::write_csv(dat, file.path(output_dir, "Fig_S05B_individual_mouse_uptake.csv"))
  readr::write_csv(sum_dat, file.path(output_dir, "Fig_S05B_mean_SEM.csv"))
  list(tables = list(uptake_summary = sum_dat, individual_uptake = dat))
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result) {
  library(ggplot2)
  dat <- result$tables$individual_uptake
  sum_dat <- result$tables$uptake_summary
  sum_dat$ymin <- pmax(0, sum_dat$mean_conc - sum_dat$se_conc)
  sum_dat$ymax <- sum_dat$mean_conc + sum_dat$se_conc
  group_levels <- levels(dat$group_label)
  base_family <- "sans"

  # open circles, faint raw dots, and SEM shown around the line.
  pal <- c(
    "IL-2/[⁸⁹Zr]hnTCB2" = "#D81B60",
    "[⁸⁹Zr]hnTCB2"      = "#1E88E5"
  )

  line_width      <- 0.58
  point_size      <- 1.08
  point_stroke    <- 0.42
  raw_point_size  <- 0.43
  raw_point_alpha <- 0.20
  sem_alpha       <- 0.13
  linerange_lwd   <- 0.34
  base_theme <- function(base_size = 6.1, base_family = "sans") {
    theme_classic(base_size = base_size, base_family = base_family) +
      theme(
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major.y = element_line(colour = "#E7E8EF", linewidth = 0.23),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.line  = element_line(colour = "#222222", linewidth = 0.33),
        axis.ticks = element_line(colour = "#222222", linewidth = 0.28),
        axis.ticks.length = unit(1.25, "pt"),
        axis.text.x = element_text(colour = "#222222", size = 5.25, margin = margin(t = 1.5)),
        axis.text.y = element_text(colour = "#222222", size = 5.25, margin = margin(r = 1.5)),
        axis.title.x = element_text(size = 6.0, margin = margin(t = 4.5)),
        axis.title.y = element_text(size = 6.0, margin = margin(r = 4.5)),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", colour = "#111111", size = 6.6, margin = margin(b = 2.5)),
        panel.spacing.x = unit(8.5, "pt"),
        panel.spacing.y = unit(7.0, "pt"),
        legend.position = "top",
        legend.justification = "right",
        legend.box.just = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 5.7, colour = "#222222"),
        legend.key = element_blank(),
        legend.key.width = unit(18, "pt"),
        legend.margin = margin(b = -2),
        plot.margin = margin(1, 4, 1, 2, "pt")
      )
  }

  facet_layer <- ggh4x::facet_wrap2(vars(organ), ncol = 4, scales = "free_y", axes = "margins")

  p_s05b <- ggplot() +
    geom_hline(yintercept = 0, colour = "#D6D8E0", linewidth = 0.25) +
    # Raw biological replicates: faint dots, jittered only slightly for visibility.
    geom_point(
      data = dat,
      aes(x = t, y = conc, colour = group_label),
      position = position_jitter(width = 0.018, height = 0, seed = 3010),
      size = raw_point_size,
      alpha = raw_point_alpha,
      stroke = 0,
      show.legend = FALSE
    )

  p_s05b <- p_s05b +
      geom_ribbon(
        data = sum_dat,
        aes(x = t, ymin = ymin, ymax = ymax, fill = group_label, group = group_label),
        alpha = sem_alpha,
        colour = NA,
        show.legend = FALSE
      )

  p_s05b <- p_s05b +
    geom_line(
      data = sum_dat,
      aes(x = t, y = mean_conc, colour = group_label, group = group_label),
      linewidth = line_width,
      lineend = "round"
    ) +
    geom_linerange(
      data = sum_dat,
      aes(x = t, ymin = ymin, ymax = ymax, colour = group_label),
      linewidth = linerange_lwd,
      alpha = 0.95,
      show.legend = FALSE
    ) +
    geom_point(
      data = sum_dat,
      aes(x = t, y = mean_conc, colour = group_label),
      shape = 21,
      fill = "white",
      size = point_size,
      stroke = point_stroke
    ) +
    facet_layer +
    scale_colour_manual(values = pal, breaks = group_levels) +
    scale_fill_manual(values = pal, breaks = group_levels) +
    scale_x_continuous(
      limits = c(-0.08, 7.12),
      breaks = 0:7,
      labels = 0:7,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = scales::breaks_pretty(n = 3),
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    guides(
      colour = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(linewidth = 0.9, size = 1.7, alpha = 1, fill = "white")
      )
    ) +
    labs(
      x = "Time (Days)",
      y = "Concentration (%ID/g)",
      caption = "Faint points: individual mice. Lines: means. Shading and bars: ± SEM. UB: urinary bladder."
    ) +
    coord_cartesian(clip = "off") +
    base_theme(base_family = base_family)
  list("Figure S5b – tissue biodistribution" = p_s05b)
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
  cat("Supplementary Fig. 5b | Tissue radiotracer uptake\nDescriptive mean and SEM across mice. Time-zero values are supplied in Source Data; the first scan was at 15 minutes.\n")
  show("uptake_summary", "Tissue uptake by tracer and time (percent injected dose/g)", c("group_label", "organ", "t", "n", "mean_conc", "sd_conc", "se_conc"))
  cat("\nComplete result tables: ", paste(names(result$tables), collapse = ", "), "\n", sep = "")
  cat("Results directory: ", result$output_dir, "\n", sep = "")
  invisible(result)
}

# Source calculates the results, prints the analysis summary and displays the plots.
if (isTRUE(getOption("SLC3010.autorun", TRUE))) {
  analysis_result <- figure_result <- NULL
  analysis_result <- local({
    needed <- c("digest", "dplyr", "readr", "readxl")
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
    paths <- c(frames[basename(frames) == "Fig_S05B_biodistribution.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_S05B_biodistribution.R"], join(getwd(), "Fig_S05B_biodistribution.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_S05B_biodistribution.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply one workbook path or directory, or use Source in RStudio.", call. = FALSE)
    path <- path.expand(if (length(args)) args[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR"))
    book <- "Source_Data.xlsx"
    sheets <- c("Supplementary Figure 5")
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
      run <- join(base, "Fig_S05B_biodistribution", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
      needed <- c("ggh4x", "ggplot2", "scales")
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
        if (is.null(size)) size <- c(11, 7)
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
