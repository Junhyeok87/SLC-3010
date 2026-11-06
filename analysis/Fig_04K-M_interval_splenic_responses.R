# Fig. 4k-m | Splenic responses by administration interval
# Fit a one-way ANOVA per readout and calculate two-sided Dunnett-adjusted comparisons of each
# regimen with the hIL-2-only reference from the individual mouse values.
# The comparison table includes both the R estimates and the reported adjusted P values used
# for the panel annotations.
#
# Input: put Source_Data.xlsx beside this script or in the repository's data/ folder.
# Numerical packages: digest, dplyr, readr, readxl, stringr, tibble.
# Plotting packages: dplyr, ggbeeswarm, ggplot2, patchwork, scales, tibble.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$comparisons
#   analysis_result$output_dir   # Folder containing Fig_04K-M_recomputed_Dunnett.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Fig. 4k-m: Interval-study splenic MP CD8 T cells, Tregs and their ratio.
# Brackets summarize the original Dunnett-adjusted comparisons against the hIL-2-only reference.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  suppressPackageStartupMessages(library(dplyr))

  fail <- function(...) stop(paste0(...), call. = FALSE)

  input_data <- readxl::read_excel(workbook, sheet = "Figure 4", range = "A781:F970", col_types = c("text", "text", "text", "numeric", "text", "numeric"), .name_repair = "minimal")

  colnames(input_data) <- stringr::str_trim(colnames(input_data))

  if (anyDuplicated(colnames(input_data))) fail("Input column names are duplicated after trimming.")

  required_cols <- c("panel_id", "metric", "group", "treatment_order", "replicate", "value")

  absent <- setdiff(required_cols, colnames(input_data))

  if (length(absent)) fail("Missing required columns: ", paste(absent, collapse = ", "))

  if (!nrow(input_data)) fail("The input table has no observations.")

  to_number <- function(x, column) {
      z <- suppressWarnings(as.numeric(as.character(x)))
      if (any(!is.finite(z))) 
          fail("Missing/non-numeric/non-finite ", column, " at input rows: ", paste(which(!is.finite(z)), collapse = ", "), 
              ". No observations were silently removed.")
      z
  }

  source_panel <- trimws(as.character(input_data$panel_id))
  panel_metrics <- c(Fig_04K = "MP_CD8T", Fig_04L = "Treg", Fig_04M = "MP_CD8T_Treg_ratio")
  if (anyNA(source_panel) || anyNA(input_data$metric) || !all(source_panel %in% names(panel_metrics)) ||
      any(as.character(input_data$metric) != unname(panel_metrics[source_panel])))
      fail("Figure 4k-m requires the matching Fig_04K, Fig_04L and Fig_04M panel IDs and readouts.")

  raw_df <- as.data.frame(input_data)

  for (column in c("panel_id", "metric", "group", "replicate")) {
      raw_df[[column]] <- trimws(as.character(raw_df[[column]]))
      if (anyNA(raw_df[[column]]) || any(!nzchar(raw_df[[column]]))) 
          fail("Empty values in column ", column, ". No observations were removed.")
  }

  raw_df$treatment_order <- to_number(raw_df$treatment_order, "treatment_order")

  raw_df$value <- to_number(raw_df$value, "value")

  if (any(raw_df$value < 0)) fail("Negative readout values found. Check the source data.")

  metric_levels <- c("MP_CD8T", "Treg", "MP_CD8T_Treg_ratio")

  experiment_levels <- c("initial_dose", "interval", "day7_dose")

  # Published panel IDs identify the experiment directly.
  infer_experiment <- function(panel_id) {
      mapping <- c(Fig_04K = "interval", Fig_04L = "interval", Fig_04M = "interval")
      unname(mapping[as.character(panel_id)])
  }

  normalize_group <- function(x) {
      x <- as.character(x)
      x <- stringr::str_trim(x)
      dplyr::case_when(stringr::str_detect(x, "^G\\d+$") ~ x, stringr::str_detect(x, "^g\\d+$") ~ stringr::str_to_upper(x), 
          stringr::str_detect(x, "^\\d+$") ~ paste0("G", x), TRUE ~ x)
  }

  # Every interval-study spleen was collected four days after hIL-2 or control administration.
  measurements <- raw_df %>% dplyr::mutate(
      experiment = factor(infer_experiment(panel_id), levels = experiment_levels),
      group = normalize_group(group), metric = factor(as.character(metric), levels = metric_levels),
      collection_day = 4L)

  if (any(is.na(measurements$experiment))) {
      print(unique(raw_df$panel_id))
      stop("Some experiment values are NA. Check panel_id naming.")
  }

  if (any(is.na(measurements$metric))) {
      print(unique(raw_df$metric))
      stop("Some metric values are not in metric_levels.")
  }

  expected_groups <- paste0("G", 1:9)

  canonical <- paste(measurements$experiment, measurements$metric, measurements$group, measurements$replicate, sep = "|")

  if (anyDuplicated(canonical)) fail("Duplicate interval metric/group/replicate keys found. ", "No duplicate biological observations were pooled or silently dropped.")

  for (met in metric_levels) {
      rows <- measurements[as.character(measurements$metric) == met, , drop = FALSE]
      if (!setequal(unique(rows$group), expected_groups)) 
          fail("Incomplete or unexpected groups for interval / ", met, ". Expected: ", paste(expected_groups, collapse = ", "))
      if (any(rows$treatment_order != match(as.character(rows$group), expected_groups))) 
          fail("treatment_order conflicts with the G1...G9 order for interval / ", met, ".")
  }

  # Original one-way ANOVA/Dunnett output. Each row is one comparison with its reference.
  dunnett <- readxl::read_excel(workbook, sheet = "Dunnett comparisons", range = "A5:F65", col_types = c(rep("text", 5), "numeric"), .name_repair = "minimal")
  if (!identical(names(dunnett), c("Figure", "Panel", "Readout", "Reference_group",
                                   "Compared_group", "Adjusted_P_value")) ||
      anyNA(dunnett) || any(!is.finite(dunnett$Adjusted_P_value)) ||
      any(dunnett$Adjusted_P_value < 0 | dunnett$Adjusted_P_value > 1) ||
      anyDuplicated(dunnett[c("Figure", "Panel", "Reference_group", "Compared_group")])) {
      stop("Check the panel, group identifiers and adjusted P values in the Dunnett comparisons worksheet.")
  }
  dunnett$Significance <- as.character(cut(dunnett$Adjusted_P_value,
      breaks = c(0, 0.0001, 0.001, 0.005, 0.05, Inf), right = FALSE,
      labels = c("****", "***", "**", "*", "ns")))

  dunnett <- dunnett[dunnett$Figure == "Figure 4", ]


  # Refit each original ordinary one-way ANOVA from the individual measurements.
  # Annotate with the reported Prism P values and export the independently
  # calculated Dunnett-adjusted P values alongside them. The treatment-minus-control effect below has an
  # explicit direction regardless of the display order of the original contrast.
  mid_labels <- c(G1 = "PBS, day 0", G2 = "hIL-2 6 \u00b5g, day 0",
      setNames(sprintf("hnTCB2 30 \u00b5g on day -%d; hIL-2 6 \u00b5g on day 0",
                       c(20, 15, 10, 5, 3, 1)), paste0("G", 3:8)),
      G9 = "hnTCB2/hIL-2 complex (30/6 \u00b5g), day 0")

  fitted <- local({
  analysis_panels <- data.frame(panel = c("k", "l", "m"), experiment = "interval", metric = metric_levels)

  anova_group_sets <- list(interval = list(control = "G2", treatments = paste0("G", c(1, 3:9))))
  # Treatment assignments and timing in Fig. 4j.
  treatment_labels <- function(experiment, groups) unname(mid_labels[groups])


  # Two-sided, single-step Dunnett probabilities. Conditional on the common
  # standard-normal control and pooled chi-squared variance, the treatment
  # normals are independent. Integrating the conditional tail directly avoids
  # losing very small P values to subtraction from one.
  quadrature_rule <- function(diagonal, off_diagonal) {
      matrix <- diag(diagonal)
      adjacent <- cbind(seq_along(off_diagonal), seq_along(off_diagonal) + 1L)
      matrix[adjacent] <- matrix[adjacent[, 2:1]] <- off_diagonal
      eig <- eigen(matrix, symmetric = TRUE)
      list(nodes = eig$values, weights = eig$vectors[1L, ]^2)
  }
  dunnett_probability <- function(abs_t, control_n, treatment_n, df, order = 64L) {
      k <- seq_len(order - 1L)
      normal <- quadrature_rule(numeric(order), sqrt(k))   # N(0,1)
      shape <- df / 2
      gamma <- quadrature_rule(2 * (seq_len(order) - 1) + shape,
                               sqrt(k * (k + shape - 1))) # Gamma(shape,1)
      control <- rep(normal$nodes / sqrt(control_n), times = order)
      scale <- rep(sqrt(pmax(0, gamma$nodes) / shape), each = order)
      weights <- as.vector(outer(normal$weights, gamma$weights))
      vapply(abs_t, function(threshold) {
          log_inside <- numeric(length(weights))
          for (n in treatment_n) {
              limit <- threshold * scale * sqrt(1 / control_n + 1 / n)
              outside <- pnorm((control - limit) * sqrt(n)) +
                         pnorm((control + limit) * sqrt(n), lower.tail = FALSE)
              log_inside <- log_inside + log1p(-pmin(1, outside))
          }
          min(1, sum(weights * (-expm1(log_inside))))
      }, numeric(1))
  }

  anova_rows <- comparison_rows <- vector("list", nrow(analysis_panels))
  for (i in seq_len(nrow(analysis_panels))) {
      panel <- analysis_panels[i, ]
      design <- anova_group_sets[[panel$experiment]]
      groups <- c(design$control, design$treatments)
      labels <- treatment_labels(panel$experiment, groups)
      if (anyNA(labels) || any(!nzchar(labels))) stop("Missing treatment labels for panel ", panel$panel)
      data <- as.data.frame(measurements[measurements$experiment == panel$experiment &
          measurements$metric == panel$metric & measurements$group %in% groups,
          c("group", "replicate", "value")])
      data$group <- factor(data$group, levels = groups)
      n <- table(data$group)
      if (any(n < 2L) || any(!is.finite(data$value)) ||
          anyNA(data$replicate) || anyDuplicated(data[c("group", "replicate")])) {
          stop("Check the individual measurements and replicates for panel ", panel$panel)
      }
      model <- lm(value ~ group, data = data)
      overall <- anova(model)
      mse <- overall[["Mean Sq"]][2]
      df <- df.residual(model)
      counts <- as.numeric(n)
      means <- as.numeric(tapply(data$value, data$group, mean))
      effects <- means[-1] - means[1]
      se <- sqrt(mse * (1 / counts[-1] + 1 / counts[1]))
      adjusted <- dunnett_probability(abs(effects / se), counts[1], counts[-1], df)
      stored <- dunnett[dunnett$Panel == panel$panel, ]
      if (nrow(stored) != length(design$treatments) ||
          !all(stored$Reference_group == design$control) ||
          !setequal(stored$Compared_group, design$treatments)) {
          stop("The reported Dunnett comparisons do not match panel ", panel$panel)
      }
      reported <- stored$Adjusted_P_value[match(design$treatments, stored$Compared_group)]
      anova_rows[[i]] <- data.frame(
          Panel = panel$panel, Metric = panel$metric, Included_groups = paste(groups, collapse = ","),
          Included_treatments = paste(labels, collapse = " | "),
          Observations = sum(counts), Model_df = overall$Df[1], Residual_df = df,
          Residual_mean_square = mse, F_statistic = overall[["F value"]][1],
          Overall_P = overall[["Pr(>F)"]][1])
      comparison_rows[[i]] <- data.frame(
          Panel = panel$panel, Metric = panel$metric, Reference_group = design$control,
          Compared_group = design$treatments, Reference_treatment = labels[1], Compared_treatment = labels[-1],
          Reference_n = counts[1], Compared_n = counts[-1],
          Reference_mean = means[1], Compared_mean = means[-1], Compared_minus_reference = effects,
          Standard_error = se, Residual_df = df, T_statistic = effects / se, Dunnett_adjusted_P = adjusted,
          Reported_adjusted_P = reported, Absolute_P_difference = abs(adjusted - reported))
  }
  anova_results <- do.call(rbind, anova_rows)
  dunnett_results <- do.call(rbind, comparison_rows)
  readr::write_csv(anova_results, file.path(output_dir, "Fig_04K-M_one_way_ANOVA.csv"))
  readr::write_csv(dunnett_results, file.path(output_dir, "Fig_04K-M_recomputed_Dunnett.csv"))
  list(anova = anova_results, comparisons = dunnett_results)
  })

  schedule <- data.frame(group = names(mid_labels), treatment = unname(mid_labels))

  summary <- measurements %>%
    group_by(experiment, metric, group) %>%
    summarise(n = n(), mean = mean(value), sd = sd(value), sem = sd / sqrt(n), .groups = "drop")
  write.csv(summary, file.path(output_dir, "Fig_04K-M_group_summary.csv"), row.names = FALSE)

  write.csv(measurements, file.path(output_dir, "Fig_04K-M_individual_measurements.csv"), row.names = FALSE)
  readr::write_csv(dunnett, file.path(output_dir, "Fig_04K-M_Dunnett_comparisons.csv"))

  # Preserve the original display: on logarithmic axes, the bars summarize
  # the log10 measurements; quantitative arithmetic means are reported separately.
  view <- measurements %>% mutate(axis_scale = if_else(as.character(metric) == "MP_CD8T_Treg_ratio", "linear", "log10"),
                           display_value = value)
  log_rows <- view$axis_scale == "log10" & view$value > 0
  view$display_value[log_rows] <- log10(view$value[log_rows])
  view$display_value[view$axis_scale == "log10" & view$value <= 0] <- NA_real_
  display_summary <- view %>% group_by(experiment, metric, group, axis_scale) %>%
    summarise(n_display = sum(is.finite(display_value)), centre = mean(display_value, na.rm = TRUE),
              spread = sd(display_value, na.rm = TRUE) / sqrt(n_display), .groups = "drop") %>%
    mutate(lower = pmax(centre - spread, 1e-06), upper = centre + spread,
           across(c(centre, lower, upper), ~ if_else(axis_scale == "log10", 10 ^ .x, .x)))

  write.csv(display_summary, file.path(output_dir, "Fig_04K-M_display_summary.csv"), row.names = FALSE)
  list(tables = list(anova = fitted$anova, comparisons = fitted$comparisons,
                     reported_comparisons = dunnett, summary = summary,
                     observations = measurements, treatment_schedule = schedule, display_summary = display_summary))
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result) {
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggbeeswarm))
  suppressPackageStartupMessages(library(scales))
  suppressPackageStartupMessages(library(grid))

  plot_data <- result$tables$observations
  display_summary <- result$tables$display_summary
  dunnett <- result$tables$reported_comparisons
  metric_levels <- levels(plot_data$metric)

  REQUESTED_FONT <- "Arial"
  fail <- function(...) stop(paste0(...), call. = FALSE)

  PLOT_FONT_FAMILY <- REQUESTED_FONT

  if (requireNamespace("systemfonts", quietly = TRUE)) {
      families <- unique(systemfonts::system_fonts()$family)
      match_font <- families[tolower(families) == tolower(REQUESTED_FONT)]
      if (length(match_font)) 
          PLOT_FONT_FAMILY <- match_font[1L]
      else {
          choices <- c("Liberation Sans", "Nimbus Sans", "DejaVu Sans")
          found <- choices[choices %in% families]
          PLOT_FONT_FAMILY <- if (length(found)) 
              found[1L]
          else "sans"
          warning(REQUESTED_FONT, " is not installed. Using ", PLOT_FONT_FAMILY, ".", call. = FALSE)
      }
  }

  metric_titles <- c(MP_CD8T = "MP CD8 T cells", Treg = "Treg", MP_CD8T_Treg_ratio = "MP CD8 T / Treg ratio")

  metric_y_labels <- list(MP_CD8T = expression("# of MP CD8 T cells (x" * 10^6 * ")"), Treg = expression("# of Tregs (x" * 
      10^6 * ")"), MP_CD8T_Treg_ratio = "MP CD8 T / Treg ratio")

  color_design <- tibble::tribble(~experiment, ~group, ~color_hex, "initial_dose", "G1", "#BDBDBD", "initial_dose", "G2", "#000000", 
      "initial_dose", "G3", "#F7A6C9", "initial_dose", "G4", "#F48AB9", "initial_dose", "G5", "#EC5AA6", "initial_dose", "G6", "#E91E63", "initial_dose", "G7", 
      "#8E1768", "initial_dose", "G8", "#214CFF", "interval", "G1", "#BDBDBD", "interval", "G2", "#000000", "interval", "G3", "#7B0F6F", 
      "interval", "G4", "#9B0F72", "interval", "G5", "#C91578", "interval", "G6", "#E91E63", "interval", "G7", "#F06FAE", 
      "interval", "G8", "#F7A6C9", "interval", "G9", "#214CFF", "day7_dose", "G1", "#BDBDBD", "day7_dose", "G2", "#E41A1C", 
      "day7_dose", "G3", "#F7A6C9", "day7_dose", "G4", "#F48AB9", "day7_dose", "G5", "#8E1768", "day7_dose", "G6", "#B9B8FF", 
      "day7_dose", "G7", "#8889FF", "day7_dose", "G8", "#5B65FF", "day7_dose", "G9", "#2748FF", "day7_dose", "G10", "#001FDB", 
      "day7_dose", "G11", "#000000")

  axis_spec <- tibble::tibble(experiment = c("initial_dose", "initial_dose", "initial_dose", "interval", "interval", "interval", "day7_dose", "day7_dose", 
  "day7_dose"), metric = c("MP_CD8T", "Treg", "MP_CD8T_Treg_ratio", "MP_CD8T", "Treg", "MP_CD8T_Treg_ratio", "MP_CD8T", 
  "Treg", "MP_CD8T_Treg_ratio"), scale_type = c("log10", "log10", "linear", "log10", "log10", "linear", "log10", 
  "log10", "linear"), y_min = c(1, 1, 0, 1, 1, 0, 1, 1, 0), y_max = c(1000, 100, 15, 150, 30, 7, 1000, 100, 20), 
      breaks = list(c(1, 10, 100, 1000), c(1, 10, 100), c(0, 5, 10, 15), c(1, 3, 10, 30, 100), c(1, 3, 10, 30), c(0, 
      1, 2, 3, 4, 5, 6, 7), c(1, 10, 100, 1000), c(1, 10, 100), c(0, 5, 10, 15, 20)))

  plot_data <- dplyr::left_join(plot_data, color_design, by = c("experiment", "group"))

  mid_labels <- setNames(result$tables$treatment_schedule$treatment, result$tables$treatment_schedule$group)

  pt_size_circle <- 3.65
  pt_stroke <- 0.75
  box_alpha <- 0.15
  box_width <- 0.52
  mean_bar_width <- 0.46
  sem_width <- 0.13
  base_font_size <- 11.2
  x_text_size_all <- 11.5
  y_text_size <- 10.4

  title_size <- 12

  theme_figure <- function(base_size = base_font_size, base_family = PLOT_FONT_FAMILY, x_text_size = x_text_size_all, 
      y_text_size_this = y_text_size, y_tick_text_size_this = y_text_size, title_size_this = title_size, x_text_angle = 0, 
      x_text_hjust = 0.5, x_text_vjust = 0.7, x_text_margin_top = 3, y_tick_margin_right = 2.2, axis_title_margin_right = 5, 
      top_margin = 8, right_margin = 12, bottom_margin = 10, left_margin = 7) theme_classic(base_size = base_size, 
      base_family = base_family) + theme(plot.title = element_text(face = "bold", size = title_size_this, hjust = 0.5, 
      margin = margin(b = 13)), plot.tag = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(color = "black", 
      size = y_text_size_this, margin = margin(r = axis_title_margin_right)), axis.text.y = element_text(color = "black", 
      size = y_tick_text_size_this, margin = margin(r = y_tick_margin_right)), axis.text.x = element_text(color = "black", 
      size = x_text_size, angle = x_text_angle, hjust = x_text_hjust, vjust = x_text_vjust, margin = margin(t = x_text_margin_top)), 
      axis.line = element_line(color = "black", linewidth = 0.4), axis.ticks = element_line(color = "black", linewidth = 0.4), 
      axis.ticks.length = unit(1.8, "mm"), legend.position = "none", panel.grid = element_blank(), plot.margin = margin(top_margin, 
          right_margin, bottom_margin, left_margin))

  get_axis_spec <- function(experiment_id, metric_id) {
      idx <- axis_spec$experiment == experiment_id & axis_spec$metric == metric_id
      if (sum(idx) != 1) 
          stop("axis_spec is missing for: ", experiment_id, " / ", metric_id)
      axis_spec[idx, , drop = FALSE]
  }

  get_y_scale <- function(experiment_id, metric_id) {
      spec <- get_axis_spec(experiment_id, metric_id)
      brks <- spec$breaks[[1]]
      if (spec$scale_type == "log10") 
          scale_y_log10(breaks = brks, labels = label_number(accuracy = 1, big.mark = ",", decimal.mark = "."), expand = expansion(mult = c(0.02, 
              0.08)))
      else scale_y_continuous(breaks = brks, labels = label_number(accuracy = 1, big.mark = ",", decimal.mark = "."), 
          expand = expansion(mult = c(0.02, 0.08)))
  }

  draw_significance_annotations <- TRUE
  sig_line_width <- 0.42
  sig_text_size <- 3.25
  sig_text_size_stars <- 3.55

  sig_ns_design <- tibble::tibble(xmin_group = "G1", xmax_group = "G2",
      y_offset_npc = 0.045, label_offset_npc = 0.018, y_line_npc = NA_real_, y_label_npc = NA_real_)
  sig_step_design <- tibble::tibble(upper_start_group = "G2", joint_group = "G7", lower_start_group = "G3",
      lower_end_group = "G9", y_upper_npc = 1.065, y_lower_npc = 1.03, label_y_npc = 1.09,
      upper_start_nudge = -0.08, joint_nudge = 0.08, lower_start_nudge = -0.08, lower_end_nudge = 0.18)

  npc_to_data_y <- function(y_npc, axis_now) {
      scale_type <- as.character(axis_now$scale_type[[1]])
      y_min <- as.numeric(axis_now$y_min[[1]])
      y_max <- as.numeric(axis_now$y_max[[1]])
      if (scale_type == "log10") {
          if (y_min <= 0 || y_max <= 0) 
              stop("Log-scale significance annotation requires positive y_min/y_max.")
          log_min <- log10(y_min)
          log_max <- log10(y_max)
          10^(log_min + y_npc * (log_max - log_min))
      }
      else y_min + y_npc * (y_max - y_min)
  }

  data_y_plus_npc <- function(y_value, offset_npc, axis_now, max_y_npc = 0.78) {
      scale_type <- as.character(axis_now$scale_type[[1]])
      y_min <- as.numeric(axis_now$y_min[[1]])
      y_max <- as.numeric(axis_now$y_max[[1]])
      if (!is.finite(y_value)) 
          return(npc_to_data_y(0.5, axis_now))
      if (scale_type == "log10") {
          if (y_value <= 0) 
              y_value <- y_min
          log_span <- log10(y_max) - log10(y_min)
          y_out <- 10^(log10(y_value) + offset_npc * log_span)
      }
      else y_out <- y_value + offset_npc * (y_max - y_min)
      pmin(y_out, npc_to_data_y(max_y_npc, axis_now))
  }

  group_to_x <- function(groups, group_levels_here, context = "significance annotation") {
      groups <- as.character(groups)
      idx <- match(groups, group_levels_here)
      if (any(is.na(idx))) {
          missing_groups <- groups[is.na(idx)]
          stop("Unknown group(s) in ", context, ": ", paste(missing_groups, collapse = ", "), "\nAvailable groups in this panel: ", 
              paste(group_levels_here, collapse = ", "))
      }
      as.numeric(idx)
  }

  add_sig_annotations <- function(p, dat_panel, metric_id, group_levels_here, axis_now, sig_text_size_this = sig_text_size, 
      sig_text_size_stars_this = sig_text_size_stars) {
      if (!isTRUE(draw_significance_annotations)) 
          return(p)
      panel <- c(MP_CD8T = "k", Treg = "l", MP_CD8T_Treg_ratio = "m")[[metric_id]]
      stats <- dunnett[dunnett$Panel == panel, ]
      if (!all(stats$Reference_group == "G2") ||
          !setequal(stats$Compared_group, setdiff(group_levels_here, "G2"))) {
          stop("Incomplete Dunnett comparisons for Figure 4", panel)
      }
      ns_now <- sig_ns_design
      ns_now$label <- stats$Significance[stats$Compared_group == "G1"]
      if (nrow(ns_now) > 0) {
          ns_anno <- ns_now %>% dplyr::rowwise() %>% dplyr::mutate(xmin = group_to_x(xmin_group, group_levels_here, 
              "sig_ns_design") - 0.18, xmax = group_to_x(xmax_group, group_levels_here, "sig_ns_design") + 0.18, 
              xmid = (xmin + xmax)/2, y_base = max(dat_panel$value[as.character(dat_panel$group) %in% c(xmin_group, 
                  xmax_group)], na.rm = TRUE), y = if (!is.na(y_line_npc)) 
                  npc_to_data_y(y_line_npc, axis_now)
              else data_y_plus_npc(y_base, offset_npc = y_offset_npc, axis_now = axis_now, max_y_npc = 0.78), y_text = if (!is.na(y_label_npc)) 
                  npc_to_data_y(y_label_npc, axis_now)
              else data_y_plus_npc(y_base, offset_npc = y_offset_npc + label_offset_npc, axis_now = axis_now, max_y_npc = 0.82)) %>% 
              dplyr::ungroup()
          p <- p + geom_segment(data = ns_anno, aes(x = xmin, xend = xmax, y = y, yend = y), inherit.aes = FALSE, 
              linewidth = sig_line_width, color = "black", lineend = "butt") + geom_text(data = ns_anno, aes(x = xmid, 
              y = y_text, label = label), inherit.aes = FALSE, size = sig_text_size_this, color = "black", family = PLOT_FONT_FAMILY, 
              vjust = 0)
      }
      step_now <- sig_step_design
      group_labels <- unique(stats$Significance[stats$Compared_group != "G1"])
      if (length(group_labels) != 1L) stop("The stepped bracket requires a common category for G3 through G9.")
      step_now$label <- group_labels
      if (nrow(step_now) > 0) {
          step_anno <- step_now %>% dplyr::rowwise() %>% dplyr::mutate(x_upper_start = group_to_x(upper_start_group, 
              group_levels_here, "sig_step_design") + upper_start_nudge, x_joint = group_to_x(joint_group, group_levels_here, 
              "sig_step_design") + joint_nudge, x_lower_start = group_to_x(lower_start_group, group_levels_here, 
              "sig_step_design") + lower_start_nudge, x_lower_end = group_to_x(lower_end_group, group_levels_here, 
              "sig_step_design") + lower_end_nudge, x_label = (x_upper_start + x_joint)/2, y_upper = npc_to_data_y(y_upper_npc, 
              axis_now), y_lower = npc_to_data_y(y_lower_npc, axis_now), y_tip = npc_to_data_y(y_lower_npc - 0.025, 
              axis_now), y_label = npc_to_data_y(label_y_npc, axis_now)) %>% dplyr::ungroup()
          step_segments <- dplyr::bind_rows(step_anno %>% dplyr::transmute(x = x_upper_start, xend = x_joint, y = y_upper, 
              yend = y_upper), step_anno %>% dplyr::transmute(x = x_joint, xend = x_joint, y = y_upper, yend = y_lower), 
              step_anno %>% dplyr::transmute(x = x_lower_start, xend = x_lower_end, y = y_lower, yend = y_lower), 
              step_anno %>% dplyr::transmute(x = x_lower_start, xend = x_lower_start, y = y_lower, yend = y_tip), 
              step_anno %>% dplyr::transmute(x = x_lower_end, xend = x_lower_end, y = y_lower, yend = y_tip))
          p <- p + geom_segment(data = step_segments, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, 
              linewidth = sig_line_width, color = "black", lineend = "butt") + geom_text(data = step_anno, aes(x = x_label, 
              y = y_label, label = label), inherit.aes = FALSE, size = sig_text_size_stars_this, color = "black", 
              family = PLOT_FONT_FAMILY, fontface = "bold")
      }
      p
  }

  make_panel <- function(dat_all, metric_id) {
      dat_all <- as.data.frame(dat_all)
      idx <- as.character(dat_all[["experiment"]]) == "interval" & as.character(dat_all[["metric"]]) == metric_id
      dat_panel <- dat_all[idx, , drop = FALSE]
      if (nrow(dat_panel) == 0) 
          stop("No data for: ", "interval", " / ", metric_id)
      dat_panel[["group"]] <- as.character(dat_panel[["group"]])
      group_map <- unique(dat_panel[, c("group", "treatment_order")])
      group_map <- group_map[order(group_map$treatment_order), , drop = FALSE]
      group_levels_here <- as.character(group_map$group)
      dat_panel[["group"]] <- factor(dat_panel[["group"]], levels = group_levels_here)
      marks <- display_summary[as.character(display_summary$metric) == metric_id & as.character(display_summary$experiment) == "interval", , drop = FALSE]
      n_by_group <- dat_panel %>% dplyr::count(group, name = "n_group")
      dat_box <- dat_panel %>% dplyr::left_join(n_by_group, by = "group") %>% dplyr::filter(n_group >= 5)
      y_lab <- metric_y_labels[[metric_id]]
      axis_now <- get_axis_spec("interval", metric_id)
      x_text_size_now <- x_text_size_all
      y_text_size_now <- y_text_size
      y_tick_size_now <- y_text_size
      title_size_now <- title_size
      sig_text_size_now <- sig_text_size
      sig_text_size_stars_now <- sig_text_size_stars
      pt_size_circle_now <- pt_size_circle
      pt_stroke_now <- pt_stroke
      right_margin_now <- 12
      left_margin_now <- 7
      axis_title_margin_right_now <- 5
      y_tick_margin_right_now <- 2.2
      x_text_margin_top_now <- 3
      x_angle_now <- 0
      x_hjust_now <- 0.5
      x_vjust_now <- 0.7
      top_margin_now <- 35
      bottom_margin_now <- 10
      p <- ggplot(dat_panel, aes(x = group, y = value)) + geom_boxplot(data = dat_box, aes(fill = color_hex, color = color_hex), 
          width = box_width, linewidth = 0.38, alpha = box_alpha, outlier.shape = NA, show.legend = FALSE) + geom_errorbar(data = marks, aes(x = group, ymin = lower, ymax = upper),
          inherit.aes = FALSE, width = sem_width, linewidth = 0.45, colour = "black") +
          geom_errorbar(data = marks, aes(x = group, ymin = centre, ymax = centre),
          inherit.aes = FALSE, width = mean_bar_width, linewidth = 0.45, colour = "black")
      p <- p + ggbeeswarm::geom_quasirandom(aes(fill = color_hex, color = color_hex), shape = 21, width = 0.18, size = pt_size_circle_now, 
          stroke = pt_stroke_now, alpha = 0.98, varwidth = FALSE, show.legend = FALSE)
      p <- add_sig_annotations(p = p, dat_panel = dat_panel, metric_id = metric_id, group_levels_here = group_levels_here, 
          axis_now = axis_now, sig_text_size_this = sig_text_size_now, sig_text_size_stars_this = sig_text_size_stars_now)
      p <- p + scale_fill_identity(guide = "none") + scale_color_identity(guide = "none") + scale_x_discrete(limits = group_levels_here, 
          drop = TRUE, expand = expansion(add = c(0.6, 0.6))) + get_y_scale("interval", metric_id) + labs(title = metric_titles[[metric_id]], 
          y = y_lab) + theme_figure(x_text_size = x_text_size_now, y_text_size_this = y_text_size_now, y_tick_text_size_this = y_tick_size_now, 
          title_size_this = title_size_now, x_text_angle = x_angle_now, x_text_hjust = x_hjust_now, x_text_vjust = x_vjust_now, 
          x_text_margin_top = x_text_margin_top_now, y_tick_margin_right = y_tick_margin_right_now, axis_title_margin_right = axis_title_margin_right_now, 
          top_margin = top_margin_now, right_margin = right_margin_now, bottom_margin = bottom_margin_now, left_margin = left_margin_now)
      if (metric_id %in% c("MP_CD8T", "Treg")) 
          p <- p + annotation_logticks(sides = "l", outside = FALSE, short = unit(0.6, "mm"), mid = unit(0.95, "mm"), 
              long = unit(1.35, "mm"), linewidth = 0.27)
      coord_y_min_now <- as.numeric(axis_now$y_min[[1]])
      coord_y_max_now <- as.numeric(axis_now$y_max[[1]])
      p <- p + coord_cartesian(ylim = c(coord_y_min_now, coord_y_max_now), clip = "off")
      p
  }

  p_mid_1 <- make_panel(dat_all = plot_data, metric_id = "MP_CD8T")

  p_mid_2 <- make_panel(dat_all = plot_data, metric_id = "Treg")

  p_mid_3 <- make_panel(dat_all = plot_data, metric_id = "MP_CD8T_Treg_ratio")

  # The nine group IDs are the original study keys. Pair them with the full dose schedule.
  key <- data.frame(group = names(mid_labels), treatment = unname(mid_labels),
                    color = plot_data$color_hex[match(names(mid_labels), as.character(plot_data$group))])
  key_plot <- ggplot(key, aes(x = 0, y = reorder(group, -seq_along(group)))) +
    geom_point(aes(colour = color), size = 2) +
    geom_text(aes(label = sprintf("%s   %s", group, treatment)), x = 0.025, hjust = 0, size = 3) +
    scale_colour_identity() + coord_cartesian(xlim = c(0, 1), clip = "off") +
    theme_void() + labs(title = "Treatment and dosing schedule") +
    theme(plot.title = element_text(size = 10, face = "bold"))
  panels <- patchwork::wrap_plots(p_mid_1, p_mid_2, p_mid_3, nrow = 1)
  combined <- patchwork::wrap_plots(panels, key_plot, ncol = 1, heights = c(2.5, 1.5)) +
    patchwork::plot_annotation(caption = "Dunnett comparisons versus G2 (hIL-2 alone).")

  list("Fig. 4k–m | Interval study and dose schedule" = combined)
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
  cat("Fig. 4k-m | Splenic responses by administration interval\nTwo-sided Dunnett comparisons with the hIL-2-only reference. Reported probabilities used in the figure and independently recomputed probabilities remain distinct.\n")
  show("anova", "One-way ANOVA", c("Panel", "Metric", "Observations", "Model_df", "Residual_df", "F_statistic", "Overall_P"))
  show("comparisons", "Dunnett comparisons: effect estimates", c("Panel", "Reference_treatment", "Compared_treatment", "Reference_n", "Compared_n", "Compared_minus_reference", "Standard_error"))
  show("comparisons", "Reported and independently recomputed adjusted P values", c("Panel", "Compared_treatment", "Dunnett_adjusted_P", "Reported_adjusted_P", "Absolute_P_difference"))
  cat("\nComplete result tables: ", paste(names(result$tables), collapse = ", "), "\n", sep = "")
  cat("Results directory: ", result$output_dir, "\n", sep = "")
  invisible(result)
}

# Source calculates the results, prints the analysis summary and displays the plots.
if (isTRUE(getOption("SLC3010.autorun", TRUE))) {
  analysis_result <- figure_result <- NULL
  analysis_result <- local({
    needed <- c("digest", "dplyr", "readr", "readxl", "stringr", "tibble")
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
    paths <- c(frames[basename(frames) == "Fig_04K-M_interval_splenic_responses.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_04K-M_interval_splenic_responses.R"], join(getwd(), "Fig_04K-M_interval_splenic_responses.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_04K-M_interval_splenic_responses.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply one workbook path or directory, or use Source in RStudio.", call. = FALSE)
    path <- path.expand(if (length(args)) args[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR"))
    book <- "Source_Data.xlsx"
    sheets <- c("Dunnett comparisons", "Figure 4")
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
      run <- join(base, "Fig_04K-M_interval_splenic_responses", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
      needed <- c("dplyr", "ggbeeswarm", "ggplot2", "patchwork", "scales", "tibble")
      missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
      if (length(missing)) stop("Install these plotting packages: ", paste(missing, collapse = ", "))
      if (interactive() && identical(Sys.getenv("RSTUDIO"), "1") && !identical(names(grDevices::dev.cur()), "RStudioGD"))
        grDevices::dev.new(noRStudioGD = FALSE)
      measure <- !interactive() && grDevices::dev.cur() == 1L
      if (measure) grDevices::pdf(NULL)
      plots <- tryCatch(suppressPackageStartupMessages(plot_results(analysis_result)), finally = if (measure) grDevices::dev.off())
      for (i in seq_along(plots)) {
        size <- analysis_result$sizes[[names(plots)[i]]]
        if (is.null(size)) size <- c(13, 6.5)
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
