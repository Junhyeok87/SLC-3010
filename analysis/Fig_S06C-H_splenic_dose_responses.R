# Supplementary Fig. 6c-h | Splenic dose-response comparisons
# Fit a one-way ANOVA for each outcome and experiment and compute two-sided, single-step
# Dunnett comparisons against the designated control.
# The lower-panel tests include the designated control and five comparison groups. The plots
# display all eleven groups.
# For panels f-h, match measurements by experiment year and original flow cytometry sample
# name. Check each recorded ratio against MP CD8 count / Treg count before the group tests.
# Both counts are reported in millions of cells; the ratio has no units.
#
# Input: put Source_Data.xlsx beside this script or in the repository's data/ folder.
# Numerical packages: digest, dplyr, readr, readxl, stringr, tibble.
# Plotting packages: dplyr, ggbeeswarm, ggplot2, patchwork, scales, stringr, tibble.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$comparisons
#   analysis_result$tables$paired_samples  # Fig_S06F-H_paired_samples.csv, one row per sample
#   analysis_result$output_dir   # Folder containing Fig_S06C-H_recomputed_Dunnett.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Supplementary Fig. 6c-h: Splenic responses to initial and day-7 hIL-2 dose titrations.
# The original Dunnett-adjusted P values for each comparison are exported alongside the figures.
# Panels f-h plot the day-4 comparator alongside the day-11 measurements.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  suppressPackageStartupMessages(library(dplyr))

  # Excel ranges used for named datasets in this analysis.
  table_ranges <- c(
    "figs06ce_initial_dose_splenic" = "'Supplementary Figure 6'!A19:F91",
    "figs06fh_day7_dose_splenic" = "'Supplementary Figure 6'!A97:I295",
    "splenic_Dunnett_comparisons" = "'Dunnett comparisons'!A5:F65"
  )

  required_cols <- c("panel_id", "metric", "group", "treatment_order", "replicate", "value")
  types <- c("text", "text", "text", "numeric", "text", "numeric")
  initial <- readxl::read_excel(workbook, range = table_ranges[["figs06ce_initial_dose_splenic"]],
                                col_types = types, .name_repair = "minimal")
  day7 <- readxl::read_excel(workbook, range = table_ranges[["figs06fh_day7_dose_splenic"]],
                             col_types = c(types, "skip", "numeric", "text"), .name_repair = "minimal")
  names(initial) <- trimws(names(initial))
  names(day7) <- trimws(names(day7))
  if (!identical(names(initial), required_cols)) {
    stop("Check the S6c-e column names in Supplementary Figure 6!A19:F91.", call. = FALSE)
  }
  if (!identical(names(day7), c(required_cols, "experiment_year", "source_sample"))) {
    stop("S6f-h requires experiment_year and source_sample in Supplementary Figure 6!H97:I295. ",
         "Select the matching Source_Data.xlsx accompanying this code release.", call. = FALSE)
  }
  panel_metrics <- c(Fig_S06C = "MP_CD8T", Fig_S06D = "Treg", Fig_S06E = "MP_CD8T_Treg_ratio",
                     Fig_S06F = "MP_CD8T", Fig_S06G = "Treg", Fig_S06H = "MP_CD8T_Treg_ratio")
  if (anyNA(initial$panel_id) || anyNA(day7$panel_id) ||
      !all(initial$panel_id %in% names(panel_metrics)[1:3]) ||
      !all(day7$panel_id %in% names(panel_metrics)[4:6]))
    stop("S6c-e and S6f-h require their matching published Fig_S06C through Fig_S06H panel IDs.", call. = FALSE)
  input_data <- dplyr::bind_rows(initial, day7)
  if (anyNA(input_data$metric) ||
      any(input_data$metric != unname(panel_metrics[input_data$panel_id])))
    stop("Each S6 panel ID must match its recorded readout.", call. = FALSE)

  raw_df <- input_data %>% dplyr::mutate(panel_id = as.character(panel_id), metric = as.character(metric), group = as.character(group), 
      treatment_order = as.numeric(treatment_order), replicate = as.numeric(replicate), value = as.numeric(value))

  metric_levels <- c("MP_CD8T", "Treg", "MP_CD8T_Treg_ratio")

  experiment_levels <- c("initial_dose", "interval", "day7_dose")

  # Published panel IDs identify the experiment directly.
  infer_experiment <- function(panel_id) {
      mapping <- c(Fig_S06C = "initial_dose", Fig_S06D = "initial_dose",
                   Fig_S06E = "initial_dose", Fig_S06F = "day7_dose",
                   Fig_S06G = "day7_dose", Fig_S06H = "day7_dose")
      unname(mapping[as.character(panel_id)])
  }

  normalize_group <- function(x) {
      x <- as.character(x)
      x <- stringr::str_trim(x)
      dplyr::case_when(stringr::str_detect(x, "^G\\d+$") ~ x, stringr::str_detect(x, "^g\\d+$") ~ stringr::str_to_upper(x), 
          stringr::str_detect(x, "^\\d+$") ~ paste0("G", x), TRUE ~ x)
  }

  day4_groups_bottom <- "G2"

  group_design <- tibble::tibble(group = paste0("G", 1:11)) %>% dplyr::mutate(analysis_day_bottom = dplyr::if_else(group %in% 
      day4_groups_bottom, "Analysis on Day4", "Analysis on Day11"))

  measurements <- raw_df %>% dplyr::mutate(experiment = infer_experiment(panel_id), group = normalize_group(group), 
      metric = as.character(metric)) %>% dplyr::left_join(group_design, by = "group") %>% dplyr::mutate(experiment = factor(experiment, levels = experiment_levels), 
      metric = factor(metric, levels = metric_levels), analysis_day = dplyr::case_when(as.character(experiment) == 
          "day7_dose" ~ analysis_day_bottom, TRUE ~ NA_character_), analysis_day = factor(analysis_day, levels = c("Analysis on Day4", 
      "Analysis on Day11"))) %>%
      dplyr::select(-analysis_day_bottom)

  if (any(is.na(measurements$experiment))) {
      stop("Some experiment values are NA. Check panel_id naming.")
  }

  if (any(is.na(measurements$metric))) {
      print(unique(raw_df$metric))
      stop("Some metric values are not in metric_levels.")
  }

  allowed <- list(initial_dose = paste0("G", 1:8), day7_dose = paste0("G", 1:11))
  for (experiment in names(allowed)) {
    rows <- measurements[measurements$experiment == experiment, ]
    for (metric in metric_levels) {
      subset <- rows[rows$metric == metric, ]
      if (!setequal(unique(subset$group), allowed[[experiment]]) ||
          anyDuplicated(subset[c("group", "replicate")]) || any(!is.finite(subset$value)))
        stop("Missing groups, invalid observations or duplicate replicate keys in ", experiment, " / ", metric)
    }
  }

  # FCS basenames recur across years. The year and sample name together identify one sample.
  # The replicate column is a display index within a treatment group.
  lower <- as.data.frame(measurements[measurements$experiment == "day7_dose", ])
  keys <- c("experiment_year", "source_sample")
  id_columns <- c("group", "replicate", keys)
  if (anyNA(lower[id_columns]) || any(!nzchar(trimws(lower$source_sample))) ||
      any(!lower$experiment_year %in% c(2023, 2024)) || any(!lower$replicate %in% 1:6)) {
    stop("Check the year, source sample, group and replicate for every S6f-h measurement.", call. = FALSE)
  }
  # The original FCS basenames start with the recorded group, e.g. G2_Day4_1.fcs.
  if (any(!grepl("^[^/\\\\]+[.]fcs$", lower$source_sample, ignore.case = TRUE)) ||
      any(sub("_.*$", "", lower$source_sample) != lower$group)) {
    stop("Check that each original S6f-h FCS basename identifies its recorded group.", call. = FALSE)
  }
  paired <- unique(lower[id_columns])
  if (anyDuplicated(lower[c(keys, "metric")]) || anyDuplicated(paired[keys]) ||
      anyDuplicated(paired[c("group", "replicate")])) {
    stop("Each S6f-h year/sample must have a single group and replicate and one value per readout.",
         call. = FALSE)
  }
  n_by_group_year <- table(factor(paired$group, levels = allowed$day7_dose),
                           factor(paired$experiment_year, levels = c(2023, 2024)))
  if (any(n_by_group_year != 3L)) {
    stop("S6f-h requires three distinct samples per treatment group in each experiment year.", call. = FALSE)
  }
  fields <- c(MP_CD8T = "mp_cd8_count", Treg = "treg_count", MP_CD8T_Treg_ratio = "recorded_ratio")
  for (metric in names(fields)) {
    values <- lower[lower$metric == metric, c(keys, "value")]
    names(values)[3L] <- fields[[metric]]
    paired <- merge(paired, values, by = keys, all = TRUE, sort = FALSE)
  }
  if (anyNA(paired[unname(fields)]) || any(paired$treg_count <= 0)) {
    stop("S6f-h requires all three measurements and a positive Treg count for each sample.",
         call. = FALSE)
  }
  paired$calculated_ratio <- paired$mp_cd8_count / paired$treg_count
  paired$ratio_difference <- paired$recorded_ratio - paired$calculated_ratio
  # The workbook reports counts and ratios at limited decimal precision.
  mismatch <- which(abs(paired$ratio_difference) > 0.00001)
  if (length(mismatch)) {
    i <- mismatch[1L]
    stop(sprintf(paste0("S6f-h ratio mismatch for %s / %s: recorded %.8g, calculated %.8g ",
                        "(absolute tolerance 0.00001). Check the source sample labels and counts."),
                 paired$experiment_year[i], paired$source_sample[i],
                 paired$recorded_ratio[i], paired$calculated_ratio[i]), call. = FALSE)
  }
  paired <- paired[order(match(paired$group, allowed$day7_dose), paired$replicate),
                   c(id_columns, unname(fields), "calculated_ratio", "ratio_difference")]
  rownames(paired) <- NULL
  readr::write_csv(paired, file.path(output_dir, "Fig_S06F-H_paired_samples.csv"))

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

  make_condition_rows <- function(experiment, condition, condition_order, groups, values) {
      if (length(groups) != length(values)) 
          stop("Condition table length mismatch for ", experiment, " / ", condition, ": length(groups) = ", length(groups), 
              ", length(values) = ", length(values))
      tibble::tibble(experiment = experiment, condition = condition, condition_order = condition_order, group = groups, 
          value = values)
  }

  top_condition_groups <- paste0("G", 1:8)

  bottom_condition_groups <- paste0("G", 1:11)

  condition_table_design <- dplyr::bind_rows(make_condition_rows(experiment = "initial_dose", condition = "hnTCB2", condition_order = 1, 
      groups = top_condition_groups, values = c("-", "+", "+", "+", "+", "+", "+", "-")), make_condition_rows(experiment = "initial_dose", 
      condition = "hIL-2", condition_order = 2, groups = top_condition_groups, values = c("-", "-", "1", "2.5", "5", 
      "10", "20", "-")), make_condition_rows(experiment = "initial_dose", condition = "SLC-3010", condition_order = 3, groups = top_condition_groups, 
      values = c("-", "-", "-", "-", "-", "-", "-", "+")), make_condition_rows(experiment = "day7_dose", condition = "SLC-3010/D0", 
      condition_order = 1, groups = bottom_condition_groups, values = c("-", "+", "+", "-", "-", "+", "+", "+", "+", 
      "+", "+")), make_condition_rows(experiment = "day7_dose", condition = "hIL-2/D7", condition_order = 2, groups = bottom_condition_groups, 
      values = c("-", "-", "-", "6", "-", "0.5", "1", "3", "6", "9", "-")), make_condition_rows(experiment = "day7_dose", 
      condition = "SLC-3010/D7", condition_order = 3, groups = bottom_condition_groups, values = c("-", "-", "-", 
      "-", "+", "-", "-", "-", "-", "-", "+")), make_condition_rows(experiment = "day7_dose", condition = "Analysis day", 
      condition_order = 4, groups = bottom_condition_groups, values = c("D11", "D4", "D11", "D11", "D11", "D11", 
      "D11", "D11", "D11", "D11", "D11")))

  # Refit each original ordinary one-way ANOVA from the individual measurements.
  # Annotate with the reported Prism P values and export the independently
  # calculated Dunnett-adjusted P values alongside them. The treatment-minus-control effect below has an
  # explicit direction regardless of the display order of the original contrast.
  fitted <- local({
  analysis_panels <- data.frame(panel = letters[3:8], experiment = rep(c("initial_dose", "day7_dose"), each = 3),
                             metric = rep(metric_levels, 2))

  anova_group_sets <- list(initial_dose = list(control = "G2", treatments = paste0("G", c(1, 3:8))),
                           day7_dose = list(control = "G3", treatments = paste0("G", 6:10)))
  # Derive each label from the same treatment assignments used beneath the plots.
  treatment_labels <- function(experiment, groups) {
      vapply(groups, function(group_id) {
          rows <- condition_table_design[condition_table_design$experiment == experiment &
              condition_table_design$group == group_id & condition_table_design$value != "-" &
              condition_table_design$condition != "Analysis day", ]
          pieces <- vapply(seq_len(nrow(rows)), function(j) {
              value <- rows$value[j]
              switch(rows$condition[j],
                  "hnTCB2" = "hnTCB2 50 \u00b5g",
                  "hIL-2" = paste0("hIL-2 ", value, " \u00b5g, 30 min after hnTCB2"),
                  "SLC-3010" = "SLC-3010",
                  "SLC-3010/D0" = "SLC-3010 72 \u00b5g, day 0",
                  "hIL-2/D7" = paste0("hIL-2 ", value, " \u00b5g, day 7"),
                  "SLC-3010/D7" = "SLC-3010 36 \u00b5g, day 7",
                  stop("Unknown treatment in the dose design: ", rows$condition[j]))
          }, character(1))
          if (length(pieces)) paste(pieces, collapse = "; ") else "Vehicle"
      }, character(1), USE.NAMES = FALSE)
  }


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
  readr::write_csv(anova_results, file.path(output_dir, "Fig_S06C-H_one_way_ANOVA.csv"))
  readr::write_csv(dunnett_results, file.path(output_dir, "Fig_S06C-H_recomputed_Dunnett.csv"))
  list(anova = anova_results, comparisons = dunnett_results)
  })

  dunnett <- dunnett[dunnett$Figure == "Supplementary Figure 6", ]
  schedule <- condition_table_design

  summary <- measurements %>%
    group_by(experiment, metric, group) %>%
    summarise(n = n(), mean = mean(value), sd = sd(value), sem = sd / sqrt(n), .groups = "drop")
  write.csv(summary, file.path(output_dir, "Fig_S06C-H_group_summary.csv"), row.names = FALSE)

  write.csv(dplyr::filter(measurements, experiment %in% c("initial_dose", "day7_dose")),
            file.path(output_dir, "Fig_S06C-H_individual_measurements.csv"), row.names = FALSE)
  readr::write_csv(dunnett[dunnett$Figure == "Supplementary Figure 6", ],
                   file.path(output_dir, "Fig_S06C-H_Dunnett_comparisons.csv"))

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

  write.csv(display_summary, file.path(output_dir, "Fig_S06C-H_display_summary.csv"), row.names = FALSE)
  list(tables = list(anova = fitted$anova, comparisons = fitted$comparisons,
                     paired_samples = paired,
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

  condition_table_design <- result$tables$treatment_schedule %>% dplyr::mutate(fill_hex = dplyr::case_when(value ==
      "-" ~ "#F7F7F7", condition == "hnTCB2" ~ "#E8EEFF", stringr::str_detect(condition, "SLC-3010") ~ "#E6E0FF", condition == 
      "hIL-2" & value == "1" ~ "#F8D8E8", condition == "hIL-2" & value == "2.5" ~ "#F4B9D4", condition == "hIL-2" & 
      value == "5" ~ "#EE8ABC", condition == "hIL-2" & value == "10" ~ "#D94F95", condition == "hIL-2" & value == 
      "20" ~ "#8E1768", condition == "hIL-2/D7" & value == "0.5" ~ "#FCE4F1", condition == "hIL-2/D7" & value == 
      "1" ~ "#F9D2E7", condition == "hIL-2/D7" & value == "3" ~ "#F2A6CE", condition == "hIL-2/D7" & value == "6" ~ 
      "#E85D9E", condition == "hIL-2/D7" & value == "9" ~ "#C73D7E", condition == "Analysis day" ~ "#F2F2F2", TRUE ~ 
      "#EFEFEF"), text_hex = dplyr::case_when(value %in% c("10", "20", "9") ~ "#FFFFFF", TRUE ~ "#111111"))

  pt_size_circle <- 3.65
  pt_size_day4_open <- 3.75
  pt_stroke <- 0.75
  pt_stroke_day4 <- 0.95
  box_alpha <- 0.15
  box_width <- 0.52
  mean_bar_width <- 0.46
  sem_width <- 0.13
  base_font_size <- 10.5
  x_text_size_all <- 10.5
  y_text_size <- 9.8
  title_size <- 11.2
  x_text_size_grid <- 14.5
  y_text_size_grid <- 13.5
  title_size_grid <- 15.6
  title_margin_bottom_grid <- 11
  condition_tile_value_size_grid <- 5
  condition_tile_analysis_value_size_grid <- 2.65
  condition_tile_row_label_size_grid <- 11

  theme_figure <- function(base_size = base_font_size, base_family = "Arial", x_text_size = x_text_size_all, y_text_size_this = y_text_size, 
      title_size_this = title_size, x_text_angle = 0, x_text_hjust = 0.5, x_text_vjust = 0.7, top_margin = 8, bottom_margin = 10, 
      hide_y_title = FALSE, title_margin_bottom = 4) theme_classic(base_size = base_size, base_family = base_family) + 
      theme(plot.title = element_text(face = "bold", size = title_size_this, hjust = 0.5, margin = margin(b = title_margin_bottom)), 
          plot.tag = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(color = if (isTRUE(hide_y_title)) "transparent" else "black", 
              size = y_text_size_this, margin = margin(r = 5)), axis.text.y = element_text(color = "black", size = y_text_size), 
          axis.text.x = element_text(color = "black", size = x_text_size, angle = x_text_angle, hjust = x_text_hjust, 
              vjust = x_text_vjust, margin = margin(t = 3)), axis.line = element_line(color = "black", linewidth = 0.4), 
          axis.ticks = element_line(color = "black", linewidth = 0.4), axis.ticks.length = unit(1.8, "mm"), legend.position = "none", 
          panel.grid = element_blank(), plot.margin = margin(top_margin, 12, bottom_margin, 7))

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
          scale_y_log10(breaks = brks, labels = label_number(accuracy = 1), expand = expansion(mult = c(0.02, 0.08)))
      else scale_y_continuous(breaks = brks, labels = label_number(accuracy = 1), expand = expansion(mult = c(0.02, 
          0.08)))
  }

  draw_significance_annotations <- TRUE
  sig_line_width <- 0.42
  sig_text_size <- 3
  sig_text_size_stars <- 3.35
  sig_text_size_grid <- 3.55
  sig_text_size_stars_grid <- 4

  sig_ns_design <- tibble::tribble(~experiment, ~metric, ~xmin_group, ~xmax_group, ~y_offset_npc, ~label_offset_npc, 
      ~y_line_npc, ~y_label_npc, "initial_dose", NA, "G1", "G2", 0.045, 0.018, NA_real_, NA_real_, "interval", NA, "G1", 
      "G2", 0.045, 0.018, NA_real_, NA_real_, "day7_dose", "MP_CD8T", "G3", "G6", NA_real_, NA_real_, 0.735, 
      0.78, "day7_dose", "Treg", "G3", "G6", NA_real_, NA_real_, 0.695, 0.745, "day7_dose", "MP_CD8T_Treg_ratio", 
      "G3", "G6", NA_real_, NA_real_, 0.7, 0.75)

  sig_step_design <- tibble::tribble(~experiment, ~metric, ~upper_start_group, ~joint_group, ~lower_start_group, 
      ~lower_end_group, ~label, ~y_upper_npc, ~y_lower_npc, ~label_y_npc, ~upper_start_nudge, ~joint_nudge, ~lower_start_nudge, 
      ~lower_end_nudge, "initial_dose", NA, "G2", "G6", "G3", "G8", "****", 0.935, 0.895, 0.96, -0.08, 0.08, -0.08, 0.06, 
      "interval", NA, "G2", "G7", "G3", "G9", "****", 1.065, 1.03, 1.09, -0.08, 0.08, -0.08, 0.18, "day7_dose", NA, "G3", 
      "G8", "G7", "G10", "****", 0.925, 0.88, 0.95, -0.08, 0.08, -0.08, 0.06)

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

  sig_rows_for_panel <- function(sig_tbl, experiment_id, metric_id) {
      dplyr::filter(sig_tbl, experiment == experiment_id, is.na(metric) | metric == metric_id)
  }

  add_sig_annotations <- function(p, dat_panel, experiment_id, metric_id, group_levels_here, axis_now, sig_text_size_this = sig_text_size, 
      sig_text_size_stars_this = sig_text_size_stars) {
      if (!isTRUE(draw_significance_annotations)) 
          return(p)
      panel <- switch(experiment_id, initial_dose = c("c", "d", "e"), interval = c("k", "l", "m"),
                      day7_dose = c("f", "g", "h"))[match(metric_id, metric_levels)]
      stats <- dunnett[dunnett$Panel == panel, ]
      if (experiment_id == "day7_dose") {
          reference <- "G3"
          compared <- paste0("G", 6:10)
          single_group <- "G6"
      } else {
          reference <- "G2"
          compared <- setdiff(group_levels_here, reference)
          single_group <- "G1"
      }
      if (!all(stats$Reference_group == reference) || !setequal(stats$Compared_group, compared)) {
          stop("Incomplete Dunnett comparisons for panel ", panel)
      }
      pair_labels <- stats$Significance[stats$Compared_group == single_group]
      if (length(pair_labels) != 1L) stop("A unique Dunnett comparison is required for the two-group bracket.")
      ns_now <- sig_rows_for_panel(sig_ns_design, experiment_id, metric_id)
      ns_now$label <- pair_labels
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
              y = y_text, label = label), inherit.aes = FALSE, size = sig_text_size_this, color = "black", family = "Arial", 
              vjust = 0)
      }
      # The published S6 panels show one stepped **** mark across groups with different
      # adjusted P values. The per-comparison CSV records each group's actual category.
      step_now <- sig_rows_for_panel(sig_step_design, experiment_id, metric_id)
      if (experiment_id == "interval") {
          shared <- unique(stats$Significance[stats$Compared_group != single_group])
          if (length(shared) != 1L) stop("The Figure 4 interval-study span requires a shared category.")
          step_now$label <- shared
      }
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
              family = "Arial", fontface = "bold")
      }
      p
  }

  make_panel <- function(dat_all, experiment_id, metric_id) {
      dat_all <- as.data.frame(dat_all)
      idx <- as.character(dat_all[["experiment"]]) == experiment_id & as.character(dat_all[["metric"]]) == metric_id
      dat_panel <- dat_all[idx, , drop = FALSE]
      if (nrow(dat_panel) == 0) 
          stop("No data for: ", experiment_id, " / ", metric_id)
      dat_panel[["group"]] <- as.character(dat_panel[["group"]])
      group_map <- unique(dat_panel[, c("group", "treatment_order")])
      group_map <- group_map[order(group_map$treatment_order), , drop = FALSE]
      group_levels_here <- as.character(group_map$group)
      dat_panel[["group"]] <- factor(dat_panel[["group"]], levels = group_levels_here)
      marks <- display_summary[as.character(display_summary$metric) == metric_id & as.character(display_summary$experiment) == experiment_id, , drop = FALSE]
      n_by_group <- dat_panel %>% dplyr::count(group, name = "n_group")
      dat_box <- dat_panel %>% dplyr::left_join(n_by_group, by = "group") %>% dplyr::filter(n_group >= 5)
      y_lab <- metric_y_labels[[metric_id]]
      axis_now <- get_axis_spec(experiment_id, metric_id)
      x_text_size_now <- x_text_size_grid
      y_text_size_now <- y_text_size_grid
      title_size_now <- title_size_grid
      sig_text_size_now <- sig_text_size_grid
      sig_text_size_stars_now <- sig_text_size_stars_grid
      if (experiment_id == "day7_dose") {
          x_angle_now <- 45
          x_hjust_now <- 1
          x_vjust_now <- 1
          top_margin_now <- 12
          bottom_margin_now <- 18
      }
      else {
          x_angle_now <- 0
          x_hjust_now <- 0.5
          x_vjust_now <- 0.7
          top_margin_now <- dplyr::if_else(experiment_id == "interval", 42, 8)
          bottom_margin_now <- 10
      }
      top_margin_now <- dplyr::case_when(experiment_id == "interval" ~ 26, experiment_id == "day7_dose" ~ 20, TRUE ~ 10)
      bottom_margin_now <- dplyr::case_when(experiment_id %in% c("initial_dose", "interval") ~ 14, experiment_id == "day7_dose" ~ 
          14, TRUE ~ 8)
      p <- ggplot(dat_panel, aes(x = group, y = value)) + geom_boxplot(data = dat_box, aes(fill = color_hex, color = color_hex), 
          width = box_width, linewidth = 0.38, alpha = box_alpha, outlier.shape = NA, show.legend = FALSE) + geom_errorbar(data = marks, aes(x = group, ymin = lower, ymax = upper),
          inherit.aes = FALSE, width = sem_width, linewidth = 0.45, colour = "black") +
          geom_errorbar(data = marks, aes(x = group, ymin = centre, ymax = centre),
          inherit.aes = FALSE, width = mean_bar_width, linewidth = 0.45, colour = "black")
      if (experiment_id == "day7_dose") {
          dat_day4 <- dat_panel[as.character(dat_panel[["analysis_day"]]) == "Analysis on Day4", , drop = FALSE]
          dat_day11 <- dat_panel[as.character(dat_panel[["analysis_day"]]) == "Analysis on Day11", , drop = FALSE]
          p <- p + ggbeeswarm::geom_quasirandom(data = dat_day4, aes(color = color_hex), shape = 21, fill = "white", 
              width = 0.18, size = pt_size_day4_open, stroke = pt_stroke_day4, alpha = 1, varwidth = FALSE, show.legend = FALSE) + 
              ggbeeswarm::geom_quasirandom(data = dat_day11, aes(fill = color_hex, color = color_hex), shape = 21, 
                  width = 0.18, size = pt_size_circle, stroke = pt_stroke, alpha = 0.98, varwidth = FALSE, show.legend = FALSE)
      }
      else p <- p + ggbeeswarm::geom_quasirandom(aes(fill = color_hex, color = color_hex), shape = 21, width = 0.18, 
          size = pt_size_circle, stroke = pt_stroke, alpha = 0.98, varwidth = FALSE, show.legend = FALSE)
      p <- add_sig_annotations(p = p, dat_panel = dat_panel, experiment_id = experiment_id, metric_id = metric_id, 
          group_levels_here = group_levels_here, axis_now = axis_now, sig_text_size_this = sig_text_size_now, sig_text_size_stars_this = sig_text_size_stars_now)
      p <- p + scale_fill_identity(guide = "none") + scale_color_identity(guide = "none") + scale_x_discrete(limits = group_levels_here, 
          drop = TRUE, expand = expansion(add = c(0.6, 0.6))) + get_y_scale(experiment_id, metric_id) + labs(title = metric_titles[[metric_id]], 
          y = y_lab) + theme_figure(x_text_size = x_text_size_now, y_text_size_this = y_text_size_now, title_size_this = title_size_now, 
          x_text_angle = x_angle_now, x_text_hjust = x_hjust_now, x_text_vjust = x_vjust_now, top_margin = top_margin_now, 
          bottom_margin = bottom_margin_now, hide_y_title = FALSE, title_margin_bottom = title_margin_bottom_grid)
      if (metric_id %in% c("MP_CD8T", "Treg")) 
          p <- p + annotation_logticks(sides = "l", outside = FALSE, short = unit(0.6, "mm"), mid = unit(0.95, "mm"), 
              long = unit(1.35, "mm"), linewidth = 0.27)
      p <- p + coord_cartesian(ylim = c(axis_now$y_min, axis_now$y_max), clip = "off")
      p
  }

  make_condition_tile_table <- function(experiment_id, group_levels_here) {
      condition_levels_here <- condition_table_design %>% dplyr::filter(experiment == experiment_id) %>% dplyr::arrange(condition_order) %>% 
          dplyr::pull(condition) %>% unique()
      if (length(condition_levels_here) == 0) 
          stop("No condition table design is defined for experiment: ", experiment_id)
      group_levels_here <- as.character(group_levels_here)
      dat_table <- condition_table_design %>% dplyr::filter(experiment == experiment_id, group %in% group_levels_here) %>% 
          dplyr::mutate(group = factor(group, levels = group_levels_here), condition = factor(condition, levels = rev(condition_levels_here)))
      missing_groups <- setdiff(group_levels_here, as.character(dat_table$group))
      if (length(missing_groups) > 0) 
          stop("Condition table is missing group(s) for ", experiment_id, ": ", paste(missing_groups, collapse = ", "))
      value_size_now <- condition_tile_value_size_grid
      row_label_size_now <- condition_tile_row_label_size_grid
      analysis_value_size_now <- condition_tile_analysis_value_size_grid
      dat_table <- dat_table %>% dplyr::mutate(value_size_mm = dplyr::if_else(as.character(condition) == "Analysis day", 
          analysis_value_size_now, value_size_now), value_fontface = dplyr::if_else(as.character(condition) == "Analysis day", 
          "plain", "bold"))
      p <- ggplot(dat_table, aes(x = group, y = condition)) + geom_tile(aes(fill = fill_hex), width = 0.88, height = 0.78, 
          color = "white", linewidth = 0.35, show.legend = FALSE) + geom_text(aes(label = value, color = text_hex, 
          size = value_size_mm, fontface = value_fontface), family = "Arial", show.legend = FALSE) + scale_fill_identity(guide = "none") + 
          scale_color_identity(guide = "none") + scale_size_identity(guide = "none") + scale_x_discrete(limits = group_levels_here, 
          drop = TRUE, expand = expansion(add = c(0.6, 0.6)), position = "top") + scale_y_discrete(drop = FALSE, labels = function(x) x) + 
          coord_cartesian(clip = "off") + labs(x = NULL, y = NULL) + theme_void(base_family = "Arial") + theme(axis.text.x = element_text(size = 9), 
          axis.text.y = element_text(color = "black", size = row_label_size_now, hjust = 1, margin = margin(r = 5), 
              lineheight = 0.92), axis.ticks = element_blank(), panel.grid = element_blank(), panel.background = element_rect(fill = "white", 
              color = NA), plot.background = element_rect(fill = "white", color = NA), plot.margin = margin(0, 12, 
              4, 12))
      p
  }

  make_triptych <- function(experiment, groups, caption) {
    panels <- lapply(metric_levels, function(metric) make_panel(plot_data, experiment, metric))
    schedule <- make_condition_tile_table(experiment, groups)
    chart_row <- patchwork::wrap_plots(panels, nrow = 1)
    patchwork::wrap_plots(chart_row, schedule, ncol = 1, heights = c(2.1, 1)) +
      patchwork::plot_annotation(caption = caption)
  }
  first <- make_triptych("initial_dose", paste0("G", 1:8),
    "hIL-2 dose in µg. Dunnett comparisons versus hnTCB2 alone (G2).")
  second <- make_triptych("day7_dose", paste0("G", 1:11),
    paste("hIL-2 dose in µg. Open points: day 4; filled points: day 11.",
          "Dunnett comparisons: G6–G10 versus G3; all 11 groups shown.", sep = "\n"))

  list("Fig. S6c–e | Initial dosing" = first, "Fig. S6f–h | Day-7 dosing" = second)
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
  cat("Supplementary Fig. 6c-h | Splenic dose-response comparisons\nTwo-sided, single-step Dunnett tests use each designated control. Reported and recomputed adjusted probabilities remain distinct; lower-panel tests use the specified six-group subset.\n")
  show("anova", "One-way ANOVA by panel", c("Panel", "Metric", "Observations", "Model_df", "Residual_df", "F_statistic", "Overall_P"))
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
    paths <- c(frames[basename(frames) == "Fig_S06C-H_splenic_dose_responses.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_S06C-H_splenic_dose_responses.R"], join(getwd(), "Fig_S06C-H_splenic_dose_responses.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_S06C-H_splenic_dose_responses.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply one workbook path or directory, or use Source in RStudio.", call. = FALSE)
    path <- path.expand(if (length(args)) args[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR"))
    book <- "Source_Data.xlsx"
    sheets <- c("Dunnett comparisons", "Supplementary Figure 6")
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
      run <- join(base, "Fig_S06C-H_splenic_dose_responses", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
      needed <- c("dplyr", "ggbeeswarm", "ggplot2", "patchwork", "scales", "stringr", "tibble")
      missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
      if (length(missing)) stop("Install these plotting packages: ", paste(missing, collapse = ", "))
      if (interactive() && identical(Sys.getenv("RSTUDIO"), "1") && !identical(names(grDevices::dev.cur()), "RStudioGD"))
        grDevices::dev.new(noRStudioGD = FALSE)
      measure <- !interactive() && grDevices::dev.cur() == 1L
      if (measure) grDevices::pdf(NULL)
      plots <- tryCatch(suppressPackageStartupMessages(plot_results(analysis_result)), finally = if (measure) grDevices::dev.off())
      for (i in seq_along(plots)) {
        size <- analysis_result$sizes[[names(plots)[i]]]
        if (is.null(size)) size <- c(12, 7)
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
