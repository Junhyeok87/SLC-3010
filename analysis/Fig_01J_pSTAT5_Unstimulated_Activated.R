# Fig. 1j | pSTAT5 responses in unstimulated and activated CD8 cells
# Fit a five-parameter logistic curve to positive-dose gMFI for each drug and activation state
# and estimate its half-response concentration in ng/mL.
# The fit_status table reports whether a fitted parameter reached a specified bound.
#
# Input: put Source_Data.xlsx beside this script or in the repository's data/ folder.
# Numerical packages: digest, dplyr, minpack.lm, readxl.
# Plotting packages: dplyr, ggplot2, scales.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$ec50
#   analysis_result$output_dir   # Folder containing Fig_01J_5PL_EC50_summary.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Figure 1j: pSTAT5 gMFI in unstimulated and activated CD8-gated samples.
# Five-parameter logistic fits use multiple starts; ribbons show +/-1.96 residual SD.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  library(dplyr)

  SUBSET_LABELS <- c(unstimulated = "Unstimulated", activated = "Activated")

  raw <- as.data.frame(readxl::read_excel(workbook, sheet = "Figure 1", range = "A44:E55", .name_repair = "minimal"))

  # The dose series and all fitted concentration values are in ng/mL.
  if (!"concentration_ng_per_mL" %in% names(raw) && "concentration_as_supplied" %in% names(raw))
      names(raw)[names(raw) == "concentration_as_supplied"] <- "concentration_ng_per_mL"
  input_columns <- c(conc = "concentration_ng_per_mL",
      fresh_Beta_IL2 = "Unstimulated_Beta_IL2_gMFI",
      `fresh_SLC-3010` = "Unstimulated_SLC3010_gMFI",
      activated_Beta_IL2 = "Activated_Beta_IL2_gMFI",
      `activated_SLC-3010` = "Activated_SLC3010_gMFI")
  missing <- setdiff(unname(input_columns), names(raw))
  if (length(missing)) stop("Missing input columns: ", paste(missing, collapse = ", "))
  raw <- raw[unname(input_columns)]
  names(raw) <- names(input_columns)
  expected <- names(input_columns)

  stopifnot(identical(names(raw), expected), nrow(raw) == 11, all(vapply(raw, is.numeric, logical(1))), all(is.finite(as.matrix(raw))), 
      all(raw$conc >= 0), !anyDuplicated(raw$conc))

  mapping <- data.frame(column = expected[-1], subset_key = c("unstimulated", "unstimulated", "activated", "activated"
  ), subset_label = c("Unstimulated", "Unstimulated", "Activated", "Activated"), drug = c("Beta IL-2", "SLC-3010", 
  "Beta IL-2", "SLC-3010"))

  plot_df <- do.call(rbind, lapply(seq_len(nrow(mapping)), function(i) {
      data.frame(dose = raw$conc, value = raw[[mapping$column[i]]], drug = mapping$drug[i], subset_key = mapping$subset_key[i], 
          subset_label = mapping$subset_label[i])
  }))

  plot_df <- plot_df[plot_df$dose > 0, ]

  plot_df$drug <- factor(plot_df$drug, levels = c("Beta IL-2", "SLC-3010"))

  plot_df$subset_key <- factor(plot_df$subset_key, levels = c("unstimulated", "activated"))

  plot_df$subset_label <- factor(plot_df$subset_label, levels = c("Unstimulated", "Activated"))

  stopifnot(nrow(plot_df) == 40)

  fit_5pl_model <- function(df) {
      df <- df %>% dplyr::filter(!is.na(dose), !is.na(value), dose > 0) %>% dplyr::arrange(dose)
      if (nrow(df) < 5) 
          stop("Not enough data points for 5PL fitting")
      y_min <- min(df$value, na.rm = TRUE)
      y_max <- max(df$value, na.rm = TRUE)
      lower_bottom <- 0
      lower_top <- max(y_max * 0.7, 1)
      upper_bottom <- y_max * 1.2
      upper_top <- y_max * 2
      start_bottom_vec <- c(y_min, max(0, y_min * 0.8), min(df$value[1], y_min))
      start_top_vec <- c(y_max, y_max * 1.05, y_max * 1.15)
      dose_min <- min(df$dose, na.rm = TRUE)
      dose_max <- max(df$dose, na.rm = TRUE)
      dose_mid <- median(df$dose, na.rm = TRUE)
      start_ec50_vec <- c(dose_mid, df$dose[which.min(abs(df$value - (y_min + y_max)/2))], sqrt(dose_min * dose_max))
      start_hill_vec <- c(0.6, 1, 1.5, 2.5)
      start_asym_vec <- c(0.5, 0.8, 1, 1.5, 2.5)
      start_grid <- expand.grid(bottom = start_bottom_vec, top = start_top_vec, ec50_param = start_ec50_vec, hill = start_hill_vec, 
          asym = start_asym_vec)
      best_model <- NULL
      best_rss <- Inf
      for (i in seq_len(nrow(start_grid))) {
          fit_i <- tryCatch(minpack.lm::nlsLM(value ~ bottom + (top - bottom)/((1 + (ec50_param/dose)^hill)^asym), 
              data = df, start = as.list(start_grid[i, ]), lower = c(bottom = lower_bottom, top = lower_top, ec50_param = dose_min/100, 
                  hill = 0.05, asym = 0.1), upper = c(bottom = upper_bottom, top = upper_top, ec50_param = dose_max * 
                  100, hill = 10, asym = 10), control = minpack.lm::nls.lm.control(maxiter = 1024)), error = function(e) NULL)
          if (!is.null(fit_i) && isTRUE(fit_i$convInfo$isConv)) {
              rss_i <- sum(stats::residuals(fit_i)^2, na.rm = TRUE)
              if (is.finite(rss_i) && rss_i < best_rss) {
                  best_rss <- rss_i
                  best_model <- fit_i
              }
          }
      }
      if (is.null(best_model)) 
          stop("5PL fitting failed for all starting values")
      best_model
  }

  # Fit each drug, subset and endpoint once. Predictions, EC50 and diagnostics
  # all use the same selected minimum-residual-sum-of-squares model.
  fit_groups <- function(data, metric) {
    keys <- expand.grid(subset_key = levels(data$subset_key), drug = levels(data$drug),
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    curves <- parameters <- diagnostics <- models <- vector("list", nrow(keys))
    names(models) <- paste(keys$drug, keys$subset_key, sep = " / ")
    for (i in seq_len(nrow(keys))) {
      drug <- keys$drug[i]
      subset <- keys$subset_key[i]
      group <- data[as.character(data$drug) == drug & as.character(data$subset_key) == subset, , drop = FALSE]
      group <- group[order(group$dose), , drop = FALSE]
      n <- nrow(group)
      group_key <- data.frame(drug = factor(drug, levels = levels(data$drug)),
          subset_key = factor(subset, levels = levels(data$subset_key)),
          subset_label = factor(SUBSET_LABELS[[subset]], levels = levels(data$subset_label)))
      model <- if (n) tryCatch(fit_5pl_model(group), error = function(e) e) else NULL
      ok <- !is.null(model) && !inherits(model, "error")
      state <- if (!n) "no observations" else if (!ok) "fit unavailable" else if
        (isTRUE(model$convInfo$isConv)) "converged" else "convergence not confirmed"
      detail <- if (!n) "" else if (!ok) conditionMessage(model) else as.character(model$convInfo$stopMessage)
      estimates <- data.frame(EC50 = NA_real_, Bottom = NA_real_, Top = NA_real_,
          Hill = NA_real_, Asymmetry = NA_real_, EC50_model_parameter = NA_real_, Half_response = NA_real_)
      at_limit <- character()
      if (ok) {
        models[[i]] <- model
        doses <- range(group$dose)
        co <- stats::coef(model)
        residual_sd <- stats::sd(stats::residuals(model), na.rm = TRUE)
        dose_grid <- data.frame(dose = 10^seq(log10(doses[1L]), log10(doses[2L]), length.out = 500L))
        fitted <- as.numeric(stats::predict(model, newdata = dose_grid))
        if (metric == "percentage") fitted <- pmin(pmax(fitted, 0), 100)
        curves[[i]] <- cbind(group_key[rep(1L, nrow(dose_grid)), , drop = FALSE], dose_grid,
            fit = fitted, lower = pmax(fitted - 1.96 * residual_sd, 0), upper = fitted + 1.96 * residual_sd)
        half <- co[["bottom"]] + (co[["top"]] - co[["bottom"]])/2
        dense <- data.frame(dose = 10^seq(log10(doses[1L]/100), log10(doses[2L]*100), length.out = 8000L))
        pred <- stats::predict(model, newdata = dense)
        if (metric == "percentage") pred <- pmin(pmax(pred, 0), 100)
        estimates[1, ] <- list(dense$dose[which.min(abs(pred - half))], co[["bottom"]], co[["top"]],
                               co[["hill"]], co[["asym"]], co[["ec50_param"]], half)
        ymax <- max(group$value)
        low <- c(bottom = 0, top = if (metric == "percentage") 1 else max(ymax * 0.7, 1),
                 ec50_param = doses[1L]/100, hill = 0.05, asym = 0.1)
        high <- c(bottom = if (metric == "percentage") 100 else ymax * 1.2,
                  top = if (metric == "percentage") 100 else ymax * 2,
                  ec50_param = doses[2L]*100, hill = 10, asym = 10)
        at_limit <- c(paste0(names(low), "_lower")[abs(co[names(low)] - low) <= 1e-8*pmax(1, abs(low))],
                      paste0(names(high), "_upper")[abs(co[names(high)] - high) <= 1e-8*pmax(1, abs(high))])
      }
      parameters[[i]] <- cbind(group_key, estimates)
      diagnostics[[i]] <- data.frame(drug, subset_key = subset, metric, observations = n,
          distinct_positive_doses = length(unique(group$dose)), fit_status = state, solver_message = detail,
          parameters_at_bound = paste(at_limit, collapse = ","))
    }
    status <- dplyr::bind_rows(diagnostics)
    status <- status[order(match(status$subset_key, levels(data$subset_key)), match(status$drug, levels(data$drug))), ]
    blank <- data[FALSE, c("drug","subset_key","subset_label","dose"), drop = FALSE]
    blank$fit <- blank$lower <- blank$upper <- numeric(0)
    list(parameters = dplyr::bind_rows(parameters), predictions = dplyr::bind_rows(blank, dplyr::bind_rows(curves)),
         status = status, models = models)
  }

  gmfi_fit <- fit_groups(plot_df, "gMFI")
  ec50 <- gmfi_fit$parameters
  ec50$EC50_analytic <- ec50$EC50_model_parameter/(2^(1/ec50$Asymmetry) - 1)^(1/ec50$Hill)
  ec50$concentration_unit <- "ng/mL"
  write.csv(ec50, file.path(output_dir, "Fig_01J_5PL_EC50_summary.csv"), row.names = FALSE)
  status <- gmfi_fit$status
  utils::write.csv(status, file.path(output_dir, "Fig_01J_5PL_fit_status.csv"), row.names = FALSE)
  list(tables = list(ec50 = ec50, fit_status = status), measurements = plot_df,
       predictions = gmfi_fit$predictions, models = gmfi_fit$models)
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result) {
  library(dplyr)
  BASE_FAMILY <- "Arial"
  BASE_SIZE <- 8.2

  DOSE_BREAKS <- c(0.001, 0.03, 0.64, 16, 400, 2000)

  DRUG_COLS <- c(`Beta IL-2` = "#2B60D6", `SLC-3010` = "#D62728")

  DRUG_SHAPES <- c(`Beta IL-2` = 16, `SLC-3010` = 15)

  GMFI_Y_UPPER <- 22000

  format_ec50_label <- function(x) ifelse(is.na(x), "EC50 = NA", paste0("EC50 = ", signif(x, 3), " ng/mL"))

  publication_theme <- function() {
      th <- ggplot2::theme_classic(base_size = BASE_SIZE, base_family = BASE_FAMILY)
      th + ggplot2::theme(plot.title = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank(), panel.border = ggplot2::element_rect(fill = NA, 
          colour = "grey25", linewidth = 0.28), axis.line = ggplot2::element_line(linewidth = 0.28, colour = "black"), 
          axis.ticks = ggplot2::element_line(linewidth = 0.25, colour = "black"), axis.ticks.length = grid::unit(1.3, 
              "pt"), axis.title = ggplot2::element_text(face = "bold", size = BASE_SIZE + 0.2, colour = "black"), 
          axis.title.x = ggplot2::element_text(face = "bold", size = BASE_SIZE + 0.2, colour = "black", margin = ggplot2::margin(t = 5)), 
          axis.text = ggplot2::element_text(size = BASE_SIZE - 0.8, colour = "black"), axis.text.x = ggplot2::element_text(size = BASE_SIZE - 
              1, margin = ggplot2::margin(t = 1.5)), strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold", 
              size = BASE_SIZE + 1.6, colour = "black", vjust = 0.5, margin = ggplot2::margin(t = 2, b = 4)), 
          strip.clip = "off", panel.spacing.x = grid::unit(2.8, "pt"), panel.grid.major = ggplot2::element_line(colour = "grey90", 
              linewidth = 0.18, linetype = "dashed"), panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom", 
          legend.text = ggplot2::element_text(size = 9),
          plot.caption = ggplot2::element_text(size = 8.5, hjust = 0, margin = ggplot2::margin(t = 5)),
          plot.margin = ggplot2::margin(t = 8, r = 6, b = 6, l = 6), plot.background = ggplot2::element_rect(fill = "white", 
              colour = NA), panel.background = ggplot2::element_rect(fill = "white", colour = NA))
  }

  make_metric_plot <- function(plot_df, pred_df, ec50_df, y_axis_label) {
      ec50_df$label <- format_ec50_label(ec50_df$EC50)
      y_upper <- GMFI_Y_UPPER
      y_scale <- ggplot2::scale_y_continuous(limits = c(-3300, y_upper), breaks = c(0, 5000, 10000, 15000, 20000), 
          labels = scales::comma, expand = ggplot2::expansion(mult = c(0.005, 0.02)))
      p <- ggplot2::ggplot()
      p <- p + ggplot2::geom_ribbon(data = pred_df, ggplot2::aes(x = dose, ymin = lower, ymax = upper, fill = drug), 
          alpha = 0.1, colour = NA, show.legend = FALSE)
      p <- p + ggplot2::geom_line(data = pred_df, ggplot2::aes(x = dose, y = fit, colour = drug), linewidth = 0.5, 
          lineend = "round") + ggplot2::geom_point(data = plot_df, ggplot2::aes(x = dose, y = value, colour = drug, 
          shape = drug), size = 1.18, stroke = 0.3, alpha = 0.96)
      ec50_label_df <- ec50_df %>% dplyr::mutate(label_y = dplyr::case_when(drug == "Beta IL-2" ~ -1450, drug == 
          "SLC-3010" ~ -2550, TRUE ~ -1450), tick_y_start = 0, tick_y_end = dplyr::case_when(drug == "Beta IL-2" ~ 
          -850, drug == "SLC-3010" ~ -1850, TRUE ~ -850))
      p <- p + ggplot2::geom_segment(data = ec50_df %>% dplyr::filter(!is.na(EC50), EC50 > 0), ggplot2::aes(x = EC50, 
          xend = EC50, y = 0, yend = y_upper * 0.97, colour = drug), linetype = "dashed", linewidth = 0.28, alpha = 0.46, 
          show.legend = FALSE) + ggplot2::geom_segment(data = ec50_label_df %>% dplyr::filter(!is.na(EC50), EC50 > 
          0), ggplot2::aes(x = EC50, xend = EC50, y = tick_y_start, yend = tick_y_end, colour = drug), linewidth = 0.25, 
          alpha = 0.75, show.legend = FALSE) + ggplot2::geom_label(data = ec50_label_df %>% dplyr::filter(!is.na(EC50), 
          EC50 > 0), ggplot2::aes(x = EC50, y = label_y, label = label, colour = drug), fill = scales::alpha("white", 
          0.92), linewidth = 0, label.padding = grid::unit(0.06, "lines"), size = 2.35, family = BASE_FAMILY, show.legend = FALSE)
      facet_layer <- ggplot2::facet_wrap(~subset_label, nrow = 1, scales = "fixed")
      p + facet_layer + ggplot2::scale_x_log10(breaks = DOSE_BREAKS, labels = DOSE_BREAKS, expand = ggplot2::expansion(mult = c(0.06, 
          0.08))) + y_scale + ggplot2::scale_colour_manual(values = DRUG_COLS, breaks = names(DRUG_COLS)) + ggplot2::scale_fill_manual(values = DRUG_COLS, 
          breaks = names(DRUG_COLS)) + ggplot2::scale_shape_manual(values = DRUG_SHAPES, breaks = names(DRUG_SHAPES)) + 
          ggplot2::labs(x = "Concentration (ng/mL)", y = y_axis_label, colour = NULL, fill = NULL, shape = NULL) + publication_theme()
  }

  p <- make_metric_plot(plot_df = result$measurements, pred_df = result$predictions, ec50_df = result$tables$ec50, y_axis_label = "pSTAT5 gMFI")
  figure <- p + ggplot2::labs(caption = "Shading: fitted curve \u00b1 1.96 residual SD (lower bound at zero).") +
      ggplot2::guides(fill = "none") + ggplot2::theme(legend.position = "bottom")
  list("Fig. 1j - pSTAT5 dose responses" = figure)
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
  cat("Fig. 1j | pSTAT5 responses in unstimulated and activated CD8 cells\nEC50 is the fitted half-response concentration; model diagnostics are reported separately.\n")
  show("ec50", "Fitted half-response concentrations (ng/mL)", c("drug", "subset_label", "EC50", "concentration_unit"))
  show("fit_status", "Model convergence and dose coverage", c("drug", "subset_key", "observations", "distinct_positive_doses", "fit_status", "parameters_at_bound"))
  cat("\nComplete result tables: ", paste(names(result$tables), collapse = ", "), "\n", sep = "")
  cat("Results directory: ", result$output_dir, "\n", sep = "")
  invisible(result)
}

# Source calculates the results, prints the analysis summary and displays the plots.
if (isTRUE(getOption("SLC3010.autorun", TRUE))) {
  analysis_result <- figure_result <- NULL
  analysis_result <- local({
    needed <- c("digest", "dplyr", "minpack.lm", "readxl")
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
    paths <- c(frames[basename(frames) == "Fig_01J_pSTAT5_Unstimulated_Activated.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_01J_pSTAT5_Unstimulated_Activated.R"], join(getwd(), "Fig_01J_pSTAT5_Unstimulated_Activated.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_01J_pSTAT5_Unstimulated_Activated.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply one workbook path or directory, or use Source in RStudio.", call. = FALSE)
    path <- path.expand(if (length(args)) args[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR"))
    book <- "Source_Data.xlsx"
    sheets <- c("Figure 1")
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
      run <- join(base, "Fig_01J_pSTAT5_Unstimulated_Activated", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
      needed <- c("dplyr", "ggplot2", "scales")
      missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
      if (length(missing)) stop("Install these plotting packages: ", paste(missing, collapse = ", "))
      if (interactive() && identical(Sys.getenv("RSTUDIO"), "1") && !identical(names(grDevices::dev.cur()), "RStudioGD"))
        grDevices::dev.new(noRStudioGD = FALSE)
      measure <- !interactive() && grDevices::dev.cur() == 1L
      if (measure) grDevices::pdf(NULL)
      plots <- tryCatch(suppressPackageStartupMessages(plot_results(analysis_result)), finally = if (measure) grDevices::dev.off())
      for (i in seq_along(plots)) {
        size <- analysis_result$sizes[[names(plots)[i]]]
        if (is.null(size)) size <- c(8, 5)
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
