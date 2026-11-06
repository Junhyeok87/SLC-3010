# Fig. 5i | Count-based C1D4 differential expression
# Calculate DESeq2 Wald statistics directly from integer RNA read counts. Apply
# the study-wide nonzero-count filter, then fit the six C1D4 samples (three mice
# per genotype). Select abs(log2 fold change) >= 1 and nominal P < 0.05.
# This script saves all fitted genes, model settings and sample size factors,
# followed by the selected-gene tables and two volcano views.
#
# Input: put Source_Data_RNA_inputs.xlsx beside this script or in the repository's data/
# folder.
# Numerical packages: DESeq2, digest, dplyr, readr, readxl, stringr.
# Plotting packages: ggplot2, ggrepel.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$differential_expression
#   analysis_result$output_dir   # Folder containing Fig_05I_results.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Figure 5i: count-based differential expression and nominal-DEG volcano plot.
# MP_CD8_raw_counts supplies the measured integer counts and gene annotations.
# GO_BP_reference and GO_reference_info supply the frozen role annotations.
# The fitted contrast is Triple hIL-2 versus wild-type mice at C1D4. The nominal
# DEG rule uses raw P values; BH-adjusted values are also computed and retained.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  # Excel ranges used for named datasets in this analysis.
  table_ranges <- c(
    "RNA_GO_reference" = "'GO_BP_reference'!A1:F7933",
    "RNA_GO_metadata" = "'GO_reference_info'!A1:B2"
  )
  read_gene_set_chunks <- function(dataset_id) {
    x <- as.data.frame(readxl::read_excel(workbook, range = table_ranges[[dataset_id]], col_types = "text", .name_repair = "minimal"), stringsAsFactors = FALSE)
    required <- c("Term", "Term_order", "Chunk", "Member_count", "Genes")
    if (!all(required %in% names(x))) stop("Incomplete membership table: ", dataset_id, call. = FALSE)
    for (n in c("Term_order", "Chunk", "Member_count")) {
      x[[n]] <- as.integer(x[[n]])
      if (anyNA(x[[n]])) stop("Invalid membership index in ", dataset_id, call. = FALSE)
    }
    if (anyDuplicated(x[c("Term_order", "Chunk")])) stop("Duplicated membership chunk.", call. = FALSE)
    x <- x[order(x$Term_order, x$Chunk), , drop = FALSE]
    members <- lapply(seq_len(nrow(x)), function(i) {
      g <- if (x$Member_count[i] == 0L) character() else strsplit(x$Genes[i], "|", fixed = TRUE)[[1]]
      if (length(g) != x$Member_count[i] || anyNA(g) || any(!nzchar(g)))
        stop("Invalid membership count in ", dataset_id, call. = FALSE)
      g
    })
    list(table = x, members = members)
  }

  read_go_reference <- function() {
    z <- read_gene_set_chunks("RNA_GO_reference")
    x <- z$table
    if (!"Source_rows" %in% names(x)) stop("Missing GO source-row indices.", call. = FALSE)
    row_ids <- lapply(seq_len(nrow(x)), function(i) {
      r <- as.integer(strsplit(x$Source_rows[i], "|", fixed = TRUE)[[1]])
      if (length(r) != x$Member_count[i] || anyNA(r)) stop("Invalid GO source-row indices.", call. = FALSE)
      r
    })
    d <- data.frame(gene_symbol = unlist(z$members, use.names = FALSE),
                    gs_name = rep(x$Term, x$Member_count),
                    source_row = unlist(row_ids, use.names = FALSE), stringsAsFactors = FALSE)
    d <- d[order(d$source_row), , drop = FALSE]
    if (!identical(d$source_row, seq_len(nrow(d)))) stop("GO source-row coverage is incomplete.", call. = FALSE)
    rownames(d) <- NULL
    d[, c("gene_symbol", "gs_name"), drop = FALSE]
  }


  FC_CUTOFF <- 1
  P_CUTOFF <- 0.05
  Y_AXIS_MAX <- 50
  FC_MIN <- 2

  RAW_P_MAX <- 0.05

  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(stringr))
  OUT_DIR <- output_dir

  flatten_lists <- function(x) {
      x <- as.data.frame(x)
      for (nm in names(x)) if (is.list(x[[nm]])) {
          x[[nm]] <- vapply(x[[nm]], function(z) paste(z, collapse = ";"), character(1))
      }
      x
  }

  write_table <- function(x, filename) {
      readr::write_csv(flatten_lists(x), file.path(OUT_DIR, filename), na = "")
  }

  get_msigdb_gobp <- function() {
      metadata <- readxl::read_excel(workbook, sheet = "GO_reference_info", range = "A1:B2", col_types = "text", .name_repair = "minimal")
      if (nrow(metadata) != 1L || !identical(metadata$source_db, "HS_C5_GO:BP_to_Mus_musculus"))
          stop("Unexpected frozen GO reference identity.")
      read_go_reference()
  }

  compute_c1d4_de <- function(workbook, output_dir, prefix = "MP_CD8_C1D4") {
    # The study preprocessing retains genes detected in every one of the 18 samples.
    # The C1D4 contrast then fits three independent mice per genotype, with size
    # factors and dispersion estimates calculated from those six samples.
    raw <- as.data.frame(readxl::read_excel(workbook, sheet = "MP_CD8_raw_counts",
                                          .name_repair = "minimal"), stringsAsFactors = FALSE)
    sample_ids <- unlist(lapply(c("Baseline", "C1D4", "C3D4"),
                               function(day) paste0(day, "_", 1:6)), use.names = FALSE)
    count_columns <- paste0(sample_ids, "_Read_Count")
    required <- c("Gene_ID", "Gene_Symbol", count_columns)
    if (anyDuplicated(names(raw)) || !all(required %in% names(raw)))
      stop("MP_CD8_raw_counts must contain Gene_ID, Gene_Symbol and the 18 named study count columns.", call. = FALSE)
    ids <- trimws(as.character(raw$Gene_ID))
    symbols <- trimws(as.character(raw$Gene_Symbol))
    if (!length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
      stop("Gene_ID must be nonmissing and unique in MP_CD8_raw_counts.", call. = FALSE)
    values <- lapply(raw[count_columns], function(x) suppressWarnings(as.numeric(as.character(x))))
    counts <- do.call(cbind, values)
    if (any(!is.finite(counts)) || any(counts < 0) || any(counts != round(counts)) ||
        any(counts > .Machine$integer.max))
      stop("RNA read counts must be finite nonnegative integers.", call. = FALSE)
    storage.mode(counts) <- "integer"
    rownames(counts) <- ids
    colnames(counts) <- sample_ids
    retained <- rowSums(counts == 0L) == 0L
    if (sum(retained) < 2L) stop("Too few genes pass the study-wide count filter.", call. = FALSE)
    selected_samples <- paste0("C1D4_", 1:6)
    model_counts <- counts[retained, selected_samples, drop = FALSE]
    samples <- data.frame(sample_id = selected_samples,
                          genotype = factor(rep(c("wild_type", "Triple_hIL2"), each = 3L),
                                            levels = c("wild_type", "Triple_hIL2")),
                          row.names = selected_samples)
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = model_counts, colData = samples,
                                        design = ~ genotype)
    dds <- DESeq2::DESeq(dds, test = "Wald", fitType = "parametric", sfType = "ratio",
                        betaPrior = FALSE, minReplicatesForReplace = 7L, useT = FALSE,
                        minmu = 0.5, parallel = FALSE)
    # alpha controls independent filtering for BH-adjusted P values. The nominal
    # DEG rule below is separately defined by raw P < 0.05 and absolute log2 FC >= 1.
    res <- as.data.frame(DESeq2::results(dds,
        contrast = c("genotype", "Triple_hIL2", "wild_type"), lfcThreshold = 0,
        altHypothesis = "greaterAbs", independentFiltering = TRUE, alpha = 0.1,
        pAdjustMethod = "BH", test = "Wald", parallel = FALSE))
    table <- data.frame(Gene_ID = rownames(res),
        Gene_Symbol = symbols[match(rownames(res), ids)],
        Mean_normalized_count = res$baseMean,
        Log2_fold_change = res$log2FoldChange,
        Signed_fold_change = ifelse(res$log2FoldChange >= 0, 2^res$log2FoldChange,
                                    -2^(-res$log2FoldChange)),
        Log2_fold_change_standard_error = res$lfcSE,
        Wald_statistic = res$stat, Raw_P_value = res$pvalue,
        BH_adjusted_P_value = res$padj, stringsAsFactors = FALSE)
    table$nominal_DEG <- is.finite(table$Raw_P_value) &
        abs(table$Log2_fold_change) >= 1 & table$Raw_P_value < 0.05
    samples$size_factor <- as.numeric(DESeq2::sizeFactors(dds))
    filter_summary <- data.frame(input_genes = nrow(counts), retained_genes = sum(retained),
                                  study_samples = ncol(counts), model_samples = ncol(model_counts),
                                  nominal_DEGs = sum(table$nominal_DEG))
    settings <- data.frame(parameter = c("design", "contrast", "prefilter", "test", "fitType",
                                         "sfType", "betaPrior", "independentFiltering", "alpha",
                                         "nominal_DEG_rule", "DESeq2_version"),
        value = c("~ genotype", "Triple_hIL2 versus wild_type at C1D4",
                  "count > 0 in every one of the 18 study samples", "Wald", "parametric",
                  "ratio", "FALSE", "TRUE", "0.1", "abs(log2FC) >= 1 and raw P < 0.05",
                  as.character(getNamespaceVersion("DESeq2"))))
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(table, file.path(output_dir, paste0(prefix, "_all_genes.csv")), row.names = FALSE, na = "")
    utils::write.csv(samples, file.path(output_dir, paste0(prefix, "_samples.csv")), row.names = FALSE)
    utils::write.csv(filter_summary, file.path(output_dir, paste0(prefix, "_filter_summary.csv")), row.names = FALSE)
    utils::write.csv(settings, file.path(output_dir, paste0(prefix, "_analysis_settings.csv")), row.names = FALSE)
    list(data = table, samples = samples, filter_summary = filter_summary, settings = settings)
  }

  select_nominal_degs <- function(source) {
      d <- source$data
      keep <- !is.na(d$pvalue_raw) & abs(d$signed_fold_change) >= FC_MIN & d$pvalue_raw < RAW_P_MAX
      d <- d[keep,,drop=FALSE]
      if(!nrow(d))stop("No genes satisfy the final nominal-DEG criteria.")
      rownames(d) <- NULL
      d
  }

  # Infer effect sizes and P values from counts; the reported DE worksheet is
  # not read. Save the complete fitted table before selecting genes for the plot.
  SOURCE_INFO <- compute_c1d4_de(workbook, output_dir, prefix = "Fig_05I")
  SOURCE_INFO$data$signed_fold_change <- SOURCE_INFO$data$Signed_fold_change
  SOURCE_INFO$data$log2FoldChange <- SOURCE_INFO$data$Log2_fold_change
  SOURCE_INFO$data$pvalue_raw <- SOURCE_INFO$data$Raw_P_value
  SOURCE_INFO$data$padj_bh <- SOURCE_INFO$data$BH_adjusted_P_value

  deg_results <- select_nominal_degs(SOURCE_INFO)

  bad <- is.na(deg_results$Gene_Symbol) | !is.finite(deg_results$log2FoldChange) | !is.finite(deg_results$pvalue_raw)

  if (any(bad)) {
      warning(sum(bad), " rows with missing/non-finite symbols, FC or P were excluded.")
  }

  deg_results <- deg_results[!bad, , drop = FALSE]

  if (!nrow(deg_results)) stop("No valid volcano rows.")


  if (any(deg_results$pvalue_raw < 0 | deg_results$pvalue_raw > 1)) stop("Raw P values must lie in [0,1].")

  positive_p <- deg_results$pvalue_raw[deg_results$pvalue_raw > 0]

  p_floor <- if (length(positive_p)) max(min(positive_p)/10, .Machine$double.xmin) else .Machine$double.xmin

  if (!length(positive_p)) warning("All raw P values are zero; inspect the upstream P-value precision.")

  deg_results <- deg_results %>% mutate(pvalue_plot = if_else(pvalue_raw <= 0, p_floor, pvalue_raw), neg_log10_p = -log10(pvalue_plot), 
      neg_log10_p_plot = neg_log10_p, is_out_of_range = neg_log10_p > Y_AXIS_MAX, is_sig_nominal = abs(log2FoldChange) >= 
          FC_CUTOFF & pvalue_raw < P_CUTOFF)

  msig_data <- get_msigdb_gobp()

  proliferation <- msig_data %>% filter(str_detect(gs_name, regex("cell_cycle|mitotic|dna_replication|chromosome_segregation|cell_division", 
      ignore_case = TRUE)))

  cytotoxicity <- msig_data %>% filter(str_detect(gs_name, regex("T_CELL_MEDIATED_CYTOTOXICITY|cytotoxic", ignore_case = TRUE)))

  effector <- msig_data %>% filter(str_detect(gs_name, regex("immune_effector|effector", ignore_case = TRUE)))

  activation <- msig_data %>% filter(str_detect(gs_name, regex("costimulation|co-stimulatory|activation", ignore_case = TRUE)))

  differentiation <- msig_data %>% filter(str_detect(gs_name, regex("differentiation|development|lineage", ignore_case = TRUE)))

  cytokine <- msig_data %>% filter(str_detect(gs_name, regex("cytokine|chemokine|interleukin|signaling|receptor signaling", 
      ignore_case = TRUE)))

  proliferation_genes <- unique(c(proliferation$gene_symbol, "Birc5", "Bub1", "Bub1b", "Ccna2", "Ccnb1", "Ccnb2", 
      "Cdk1", "Mki67", "Top2a", "Aurkb", "Aurka", "Plk1", "Nusap1", "Cenpf", "Ube2c", "Nuf2", "Hmmr"))

  effector_activation_genes <- unique(c(cytotoxicity$gene_symbol, effector$gene_symbol, activation$gene_symbol, "Gzma", 
      "Gzmb", "Gzmk", "Cd28", "Icos", "Cd40lg", "Tnfrsf9", "Tnfrsf4", "Lag3", "Pvrig"))

  differentiation_genes <- unique(c(differentiation$gene_symbol, "Tcf7", "Il7r", "Bcl6", "Id3", "Zbtb16", "Slamf6", 
      "Nr4a2"))

  cytokine_genes <- unique(c(cytokine$gene_symbol, "Ifng", "Il2", "Tnf", "Cxcl9", "Cxcl10", "Cxcl11", "Il12rb2"))

  assign_role <- function(gene) case_when(gene %in% proliferation_genes ~ "Proliferation", gene %in% effector_activation_genes ~ 
      "Effector/Activation", gene %in% differentiation_genes ~ "Differentiation", gene %in% cytokine_genes ~ "Cytokine/Signaling", 
      TRUE ~ "Others")

  role_levels <- c("Proliferation", "Effector/Activation", "Differentiation", "Cytokine/Signaling", "Others")

  deg_results <- deg_results %>% mutate(role = factor(assign_role(Gene_Symbol), levels = role_levels))
  # The article labels eight genes. Use data coordinates so every selected point is visible at any size.
  label_genes <- c("Hmmr", "Nuf2", "Mki67", "Aurkb", "Cdk1", "Ccnb2", "Birc5", "Top2a")
  deg_results$label <- ifelse(deg_results$Gene_Symbol %in% label_genes, deg_results$Gene_Symbol, "")
  write_table(deg_results[, c("Gene_ID", "Gene_Symbol", "signed_fold_change", "log2FoldChange", "pvalue_raw",
      "padj_bh", "neg_log10_p", "neg_log10_p_plot", "role", "label", "is_out_of_range")], "Fig_05I_results.csv")
  # Roles can overlap; the plotted role uses the original priority, in column order.
  roles <- data.frame(
      Gene_ID = deg_results$Gene_ID, Gene_Symbol = deg_results$Gene_Symbol,
      proliferation = deg_results$Gene_Symbol %in% proliferation_genes,
      effector_activation = deg_results$Gene_Symbol %in% effector_activation_genes,
      differentiation = deg_results$Gene_Symbol %in% differentiation_genes,
      cytokine_signaling = deg_results$Gene_Symbol %in% cytokine_genes,
      displayed_role = as.character(deg_results$role),
      label_selected = nzchar(deg_results$label))
  roles$n_roles <- rowSums(roles[, c("proliferation", "effector_activation", "differentiation", "cytokine_signaling")])
  write_table(roles, "Fig_05I_role_membership.csv")
  summary <- data.frame(source_rows = nrow(SOURCE_INFO$data), selected_rows = nrow(deg_results),
                        selected_bh_below_0_05 = sum(deg_results$padj_bh < 0.05, na.rm = TRUE),
                        selected_bh_unavailable = sum(is.na(deg_results$padj_bh)),
                        overlapping_roles = sum(roles$n_roles > 1L),
                        above_article_y_limit = sum(deg_results$is_out_of_range))
  write_table(summary, "Fig_05I_selection_summary.csv")
  list(tables = list(selection_summary = summary, differential_expression = deg_results,
                     all_fitted_genes = SOURCE_INFO$data, model_samples = SOURCE_INFO$samples,
                     filtering = SOURCE_INFO$filter_summary, model_settings = SOURCE_INFO$settings,
                     role_membership = roles),
       threshold = c(abs_log2_fold_change = FC_CUTOFF, raw_p = P_CUTOFF, article_y_limit = Y_AXIS_MAX))
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result) {
  suppressPackageStartupMessages(library(ggplot2))
  deg_results <- result$tables$differential_expression
  FC_CUTOFF <- unname(result$threshold["abs_log2_fold_change"])
  P_CUTOFF <- unname(result$threshold["raw_p"])
  Y_AXIS_MAX <- unname(result$threshold["article_y_limit"])
  role_levels <- levels(deg_results$role)
  role_colors <- c(Proliferation = "#E64B35", `Effector/Activation` = "#4DBBD5", Differentiation = "#8DAA3B", `Cytokine/Signaling` = "#D7A83B", 
      Others = "#B8B8B8")

  role_shapes <- c(Proliferation = 15, `Effector/Activation` = 18, Differentiation = 16, `Cytokine/Signaling` = 17, 
      Others = 16)

  role_sizes <- c(Proliferation = 1.55, `Effector/Activation` = 1.5, Differentiation = 1.45, `Cytokine/Signaling` = 1.5, 
      Others = 0.5)
  label_rows <- deg_results[nzchar(deg_results$label), , drop = FALSE]
  full_range_count <- sum(deg_results$is_out_of_range)
  caption <- sprintf(paste0("Signed fold change of at least 2 in either direction and nominal P < 0.05: ",
                             "%s of %s rows.
  Full range includes %s above the article axis limit of 50."),
                     nrow(deg_results), result$tables$selection_summary$source_rows, full_range_count)
  volcano <- ggplot(deg_results, aes(log2FoldChange, neg_log10_p_plot)) +
      geom_hline(yintercept = -log10(P_CUTOFF), linewidth = 0.3, linetype = "22", colour = "grey50") +
      geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF), linewidth = 0.3, linetype = "22", colour = "grey50") +
      geom_point(data = subset(deg_results, role == "Others"), aes(colour = role, shape = role, size = role),
                 alpha = 0.55) +
      geom_point(data = subset(deg_results, role != "Others"), aes(colour = role, shape = role, size = role),
                 alpha = 0.88) +
      ggrepel::geom_text_repel(data = label_rows, aes(label = label), size = 3, fontface = "italic",
          seed = 3010, min.segment.length = 0, max.overlaps = Inf, box.padding = 0.35, show.legend = FALSE) +
      scale_colour_manual("Gene role", values = role_colors, breaks = role_levels, drop = FALSE) +
      scale_shape_manual("Gene role", values = role_shapes, breaks = role_levels, drop = FALSE) +
      scale_size_manual("Gene role", values = role_sizes * 1.4, breaks = role_levels, drop = FALSE) +
      scale_y_continuous(expand = expansion(mult = c(0.03, 0.10))) +
      labs(x = expression(log[2]~"fold change (Triple hIL-2 / WT at C1D4)"),
           y = expression(-log[10]~"nominal P"), caption = caption) +
      guides(colour = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE),
             size = guide_legend(nrow = 2, byrow = TRUE)) +
      theme_classic(base_size = 10) + theme(legend.position = "bottom", plot.caption = element_text(hjust = 0))
  overflow <- deg_results[deg_results$is_out_of_range, , drop = FALSE]
  detail <- volcano + coord_cartesian(ylim = c(0, Y_AXIS_MAX)) +
      labs(caption = sprintf("Nominal DEG subset: %s of %s rows. Open triangles mark values above %s; see the full-range plot.",
                             nrow(deg_results), result$tables$selection_summary$source_rows, Y_AXIS_MAX))
  if (nrow(overflow)) {
      detail <- detail + geom_point(data = overflow, aes(x = log2FoldChange), y = Y_AXIS_MAX,
          inherit.aes = FALSE, shape = 24, colour = "black", fill = "white", size = 3, show.legend = FALSE) +
          ggrepel::geom_text_repel(data = overflow, aes(x = log2FoldChange, label = Gene_Symbol),
              y = Y_AXIS_MAX, inherit.aes = FALSE, seed = 3010, size = 3, nudge_y = -2, show.legend = FALSE)
  }
  list("Figure 5i | MP-CD8 volcano at article scale" = detail,
       "Figure 5i | MP-CD8 volcano at full range" = volcano)
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
  cat("Fig. 5i | Count-based C1D4 differential expression\nDESeq2 Wald test on raw counts; displayed selection uses absolute log2 fold change >= 1 and nominal P < 0.05. Full gene-level estimates and adjusted P values are saved in the result tables.\n")
  show("filtering", "Count filtering and model sample sizes", c("input_genes", "retained_genes", "study_samples", "model_samples", "nominal_DEGs"))
  show("selection_summary", "Gene selection summary", c("source_rows", "selected_rows", "selected_bh_below_0_05", "selected_bh_unavailable", "overlapping_roles", "above_article_y_limit"))
  show("model_samples", "DESeq2 sample size factors", c("sample_id", "genotype", "size_factor"))
  cat("\nComplete result tables: ", paste(names(result$tables), collapse = ", "), "\n", sep = "")
  cat("Results directory: ", result$output_dir, "\n", sep = "")
  invisible(result)
}

# Source calculates the results, prints the analysis summary and displays the plots.
if (isTRUE(getOption("SLC3010.autorun", TRUE))) {
  analysis_result <- figure_result <- NULL
  analysis_result <- local({
    needed <- c("DESeq2", "digest", "dplyr", "readr", "readxl", "stringr")
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
    paths <- c(frames[basename(frames) == "Fig_05I_volcano.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_05I_volcano.R"], join(getwd(), "Fig_05I_volcano.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_05I_volcano.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply one workbook path or directory, or use Source in RStudio.", call. = FALSE)
    path <- path.expand(if (length(args)) args[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR"))
    book <- "Source_Data_RNA_inputs.xlsx"
    sheets <- c("GO_BP_reference", "GO_reference_info", "MP_CD8_raw_counts")
    expected <- "2779beb9064ccd6090031bb707590c4de1d32a53e7b71326388a806ec7e9fd26"
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
      run <- join(base, "Fig_05I_volcano", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
      needed <- c("ggplot2", "ggrepel")
      missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
      if (length(missing)) stop("Install these plotting packages: ", paste(missing, collapse = ", "))
      if (interactive() && identical(Sys.getenv("RSTUDIO"), "1") && !identical(names(grDevices::dev.cur()), "RStudioGD"))
        grDevices::dev.new(noRStudioGD = FALSE)
      measure <- !interactive() && grDevices::dev.cur() == 1L
      if (measure) grDevices::pdf(NULL)
      plots <- tryCatch(suppressPackageStartupMessages(plot_results(analysis_result)), finally = if (measure) grDevices::dev.off())
      for (i in seq_along(plots)) {
        size <- analysis_result$sizes[[names(plots)[i]]]
        if (is.null(size)) size <- c(9, 6)
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
