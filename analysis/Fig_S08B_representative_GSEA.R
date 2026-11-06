# Supplementary Fig. 8b | Representative enrichment curves
# Refit the C1D4 comparison from raw counts, map selected genes using the supplied
# reference, and calculate pathway enrichment with frozen GO memberships.
# The figure's NES and P values are calculated in this run; no saved statistics supply them.
# Input: Source_Data_RNA_inputs.xlsx in data/ or beside this script.
# Numerical packages: readxl, digest, DESeq2, fgsea, BiocParallel.
# Plotting packages: ggplot2, patchwork, svglite.
# SVG output: 171.979166666667 x 89.1645833333333 mm; Arial text.
# Click Source in RStudio for numerical results followed by plots with legends.
# Inspect analysis_result$tables$representative_pathways and analysis_result$output_dir.
# options(SLC3010.tables_only = TRUE) skips plots; reset it to NULL for plots.
# Genes passing abs(log2FC) >= 1 and nominal P < 0.05 form the signed-linear FC ranking.
# Frozen identifier mappings and S8_pathway_annotations provide reference annotations only.
# S8_GSEA_results and S8_GSEA_all_tested are reported outputs, not analysis inputs.
# A fixed random seed is used for the enrichment calculation.

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


read_frozen_gene_sets <- function(workbook, sheet) {
  x <- as.data.frame(readxl::read_excel(workbook, sheet = sheet, col_types = "text",
                                      .name_repair = "minimal"), stringsAsFactors = FALSE)
  required <- c("Term", "Term_order", "Chunk", "Member_count", "Genes")
  if (!all(required %in% names(x))) stop("Incomplete frozen gene-set table: ", sheet, call. = FALSE)
  for (n in c("Term_order", "Chunk", "Member_count")) x[[n]] <- suppressWarnings(as.integer(x[[n]]))
  if (anyNA(x[required[2:4]]) || any(x$Member_count < 0L) ||
      anyNA(x$Term) || any(!nzchar(x$Term)) || anyDuplicated(x[c("Term_order", "Chunk")]))
    stop("Invalid frozen gene-set chunk index.", call. = FALSE)
  x <- x[order(x$Term_order, x$Chunk), , drop = FALSE]
  chunks <- lapply(seq_len(nrow(x)), function(i) {
    genes <- if (x$Member_count[i] == 0L) character() else strsplit(x$Genes[i], "|", fixed = TRUE)[[1]]
    if (length(genes) != x$Member_count[i] || anyNA(genes) || any(!nzchar(genes)))
      stop("Invalid frozen gene-set membership count.", call. = FALSE)
    genes
  })
  terms <- unique(x$Term_order)
  sets <- lapply(terms, function(term) {
    ix <- which(x$Term_order == term)
    if (length(unique(x$Term[ix])) != 1L || !identical(x$Chunk[ix], seq_along(ix)))
      stop("Inconsistent frozen gene-set chunk ordering.", call. = FALSE)
    unique(unlist(chunks[ix], use.names = FALSE))
  })
  names(sets) <- vapply(terms, function(term) x$Term[match(term, x$Term_order)], "")
  if (anyDuplicated(names(sets))) stop("Duplicated gene-set name.", call. = FALSE)
  sets
}

calculate_enrichment <- function(pathways, ranks, seed, eps, exponent, min_size, max_size, output_dir, prefix) {
  if (!length(ranks) || any(!is.finite(ranks)) || anyNA(names(ranks)) ||
      any(!nzchar(names(ranks))) || anyDuplicated(names(ranks)) || any(diff(ranks) > 0))
    stop("The enrichment input must be a finite decreasing vector with unique gene identifiers.", call. = FALSE)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  previous_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  previous_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(previous_kind))
    if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv) else
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  utils::write.csv(data.frame(index = seq_along(.Random.seed), value = .Random.seed),
                   file.path(output_dir, paste0(prefix, "_initial_rng_state.csv")), row.names = FALSE)
  result <- as.data.frame(fgsea::fgseaMultilevel(pathways = pathways, stats = ranks,
      minSize = min_size, maxSize = max_size, eps = eps, gseaParam = exponent,
      sampleSize = 101, nPermSimple = 1000, scoreType = "std",
      BPPARAM = BiocParallel::SerialParam()), stringsAsFactors = FALSE)
  result$padj <- stats::p.adjust(result$pval, method = "BH")
  settings <- data.frame(parameter = c("seed", "RNG", "parallel_backend", "minSize", "maxSize",
      "eps", "exponent", "sampleSize", "nPermSimple", "scoreType", "pAdjustMethod", "fgsea_version"),
      value = c(seed, "Mersenne-Twister / Inversion / Rejection", "BiocParallel::SerialParam",
                min_size, max_size, eps, exponent, 101, 1000, "std", "BH",
                as.character(getNamespaceVersion("fgsea"))))
  utils::write.csv(settings, file.path(output_dir, paste0(prefix, "_enrichment_settings.csv")), row.names = FALSE)
  writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, paste0(prefix, "_session_info.txt")))
  list(data = result, settings = settings)
}

write_result_table <- function(x, output_dir, name) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  for (column in names(x)) if (is.list(x[[column]]))
    x[[column]] <- vapply(x[[column]], paste, "", collapse = ";")
  utils::write.csv(x, file.path(output_dir, name), row.names = FALSE, na = "")
  invisible(x)
}


analyze_data <- function(workbook, output_dir) {
  de <- compute_c1d4_de(workbook, output_dir, "Fig_S08B_C1D4_DE")
  selected_genes <- de$data[de$data$nominal_DEG, , drop = FALSE]
  if (!nrow(selected_genes) || anyNA(selected_genes$Gene_Symbol) || any(!nzchar(selected_genes$Gene_Symbol)))
    stop("Selected differential-expression rows require gene symbols.", call. = FALSE)
  mapping <- as.data.frame(readxl::read_excel(workbook, sheet = "S8_gene_ID_mapping", range = "A1:B416",
      col_types = "text", .name_repair = "minimal"), stringsAsFactors = FALSE)
  if (!all(c("Gene_Symbol", "Entrez_ID") %in% names(mapping)) || anyNA(mapping) ||
      any(!nzchar(mapping$Gene_Symbol)) || any(!nzchar(mapping$Entrez_ID)) ||
      anyDuplicated(mapping[c("Gene_Symbol", "Entrez_ID")])) stop("Invalid identifier mapping.", call. = FALSE)
  unmapped <- selected_genes[!selected_genes$Gene_Symbol %in% mapping$Gene_Symbol, , drop = FALSE]
  input <- data.frame(Gene_Symbol = selected_genes$Gene_Symbol, ranking_metric = selected_genes$Signed_fold_change)
  mapped <- merge(mapping, input, by = "Gene_Symbol")
  ranks <- sort(setNames(mapped$ranking_metric, mapped$Entrez_ID), decreasing = TRUE)
  if (anyDuplicated(names(ranks))) stop("Mapped genes contain duplicated Entrez identifiers.", call. = FALSE)
  sets <- read_frozen_gene_sets(workbook, "S8_GSEA_gene_sets")
  pars <- as.data.frame(readxl::read_excel(workbook, sheet = "S8_GSEA_parameters", range = "A1:F2",
      .name_repair = "minimal"))
  required <- c("pvalueCutoff", "eps", "pAdjustMethod", "exponent", "minGSSize", "maxGSSize")
  if (nrow(pars) != 1L || !all(required %in% names(pars)) || pars$pAdjustMethod != "BH")
    stop("Invalid frozen enrichment parameters.", call. = FALSE)
  fit <- calculate_enrichment(sets, ranks, seed = 20261102, eps = pars$eps, exponent = pars$exponent,
      min_size = pars$minGSSize, max_size = pars$maxGSSize, output_dir, "Fig_S08B")
  calculated <- fit$data
  # Read all reference labels independently of reported results; statistics are calculated above.
  labels <- as.data.frame(readxl::read_excel(workbook, sheet = "S8_pathway_annotations",
      col_types = "text", .name_repair = "minimal"), stringsAsFactors = FALSE)
  if (!all(c("ID", "Description") %in% names(labels)) || anyDuplicated(labels$ID))
    stop("Invalid pathway-label reference.", call. = FALSE)
  all_results <- data.frame(ID = calculated$pathway,
      Description = labels$Description[match(calculated$pathway, labels$ID)],
      setSize = calculated$size, enrichmentScore = calculated$ES,
      NES = calculated$NES, pvalue = calculated$pval, p.adjust = calculated$padj,
      log2err = calculated$log2err,
      core_enrichment = vapply(calculated$leadingEdge, paste, "", collapse = "/"))
  no_label <- is.na(all_results$Description) | !nzchar(all_results$Description)
  all_results$Description[no_label] <- all_results$ID[no_label]
  # gseGO's reported table applies pvalueCutoff to both nominal and adjusted P values.
  eligible <- is.finite(all_results$pvalue) & is.finite(all_results$p.adjust) &
      all_results$pvalue <= pars$pvalueCutoff & all_results$p.adjust <= pars$pvalueCutoff
  r <- all_results[eligible, , drop = FALSE]
  r <- r[order(-r$NES, r$ID), , drop = FALSE]
  ranked <- data.frame(rank = seq_along(ranks), gene_id = names(ranks), score = unname(ranks))
  # Save numerical results before the representative-panel selection, so that an
  # unsupported selection in a different dataset cannot discard the analysis.
  write_result_table(all_results, output_dir, "Fig_S08B_all_tested_pathways.csv")
  write_result_table(r, output_dir, "Fig_S08B_significant_pathways.csv")
  write_result_table(ranked, output_dir, "Fig_S08B_ranked_genes.csv")
  write_result_table(unmapped, output_dir, "Fig_S08B_unmapped_selected_genes.csv")
  anchor <- "GO:0051301"
  top <- head(r$ID[r$NES > 0], 20L)
  if (!anchor %in% top || length(setdiff(top, anchor)) < 2L)
    stop("The calculated results do not support the specified representative-pathway selection.", call. = FALSE)
  le <- setNames(strsplit(r$core_enrichment, "/", fixed = TRUE), r$ID)
  jaccard <- function(a, b) length(intersect(le[[a]], le[[b]])) / length(union(le[[a]], le[[b]]))
  pairs <- combn(setdiff(top, anchor), 2L)
  overlap <- apply(pairs, 2L, function(pair) mean(c(jaccard(anchor, pair[1]),
      jaccard(anchor, pair[2]), jaccard(pair[1], pair[2]))))
  ids <- c(anchor, pairs[, which.min(overlap)])
  representative <- r[match(ids, r$ID), , drop = FALSE]
  curves <- lapply(seq_along(ids), function(i) {
    hit <- names(ranks) %in% sets[[ids[i]]]
    weight <- abs(ranks)^pars$exponent
    if (!any(hit) || all(hit) || sum(weight[hit]) == 0) stop("Degenerate running-score input.", call. = FALSE)
    curve <- cumsum(ifelse(hit, weight / sum(weight[hit]), -1 / sum(!hit)))
    peak <- which.max(curve)
    leading <- names(ranks)[hit & seq_along(ranks) <= peak]
    if (sum(hit) != representative$setSize[i] ||
        abs(curve[peak] - representative$enrichmentScore[i]) > 1e-7 ||
        !setequal(leading, le[[ids[i]]])) stop("Running score disagrees with calculated enrichment.", call. = FALSE)
    data.frame(pathway = ids[i], rank = c(0, seq_along(ranks)), ES = c(0, curve),
               set_member = c(FALSE, hit))
  })
  curves <- do.call(rbind, curves)
  write_result_table(representative, output_dir, "Fig_S08B_representative_pathways.csv")
  write_result_table(curves, output_dir, "Fig_S08B_running_scores.csv")
  list(tables = list(representative_pathways = representative, all_tested_pathways = all_results,
      significant_pathways = r, running_scores = curves, differential_expression = de$data,
      unmapped_selected_genes = unmapped), ranked_genes = ranked, settings = fit$settings)
}

plot_results <- function(result) {
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(patchwork))
  selected <- result$tables$representative_pathways
  curves <- result$tables$running_scores
  ranks <- result$ranked_genes$score
  ids <- selected$ID
  fmt <- function(x) formatC(x, format = "e", digits = 2)

  theme_base <- theme_classic(base_size = 7, base_family = "Arial") + theme(panel.border = element_rect(fill = NA, 
      colour = "grey45", linewidth = 0.35), axis.line = element_blank(), axis.text = element_text(size = 6.5, 
      colour = "black"), plot.margin = margin(2, 3, 2, 3), axis.title = element_text(size = 7))

  plots <- list()

  for (i in seq_along(ids)) {
      id <- ids[i]
      row <- selected[i, ]

      d <- curves[curves$pathway == id, , drop = FALSE]
      hit_ranks <- d$rank[d$set_member]
      title <- tools::toTitleCase(row$Description)
      title <- paste(strwrap(title, width = 27), collapse = "\n")
      p_floor <- as.numeric(result$settings$value[result$settings$parameter == "eps"])
      p_relation <- if (is.finite(row$pvalue) && row$pvalue <= p_floor) "\u2264" else "="
      ann <- sprintf("NES = %.3f\nNominal P %s %s\nAdjusted P = %s",
                     row$NES, p_relation, fmt(row$pvalue), fmt(row$p.adjust))
      upper <- ggplot(d, aes(rank, ES)) + geom_line(colour = "#D94C43", linewidth = 0.4) + geom_hline(yintercept = 0, 
          linetype = "dashed", linewidth = 0.3) + annotate("text", x = length(ranks) * 0.04, y = 0.235, 
          label = ann, hjust = 0, vjust = 1, size = 7/ggplot2::.pt, family = "Arial") + scale_x_continuous(limits = c(0, 
          length(ranks)), expand = c(0, 0)) + scale_y_continuous(limits = c(-0.025, 0.8), breaks = c(0, 
          0.2, 0.4, 0.6, 0.8), expand = c(0, 0)) + labs(title = title, y = "Running enrichment score", 
          x = NULL) + theme_base + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          plot.title = element_text(size = 11, face = "bold", hjust = 0.5))
      ticks <- ggplot(data.frame(rank = hit_ranks), aes(x = rank)) + geom_segment(aes(xend = rank, y = 0, 
          yend = 1), linewidth = 0.23) + scale_x_continuous(limits = c(0, length(ranks)), expand = c(0, 
          0)) + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + labs(x = NULL, y = NULL) + theme_base + 
          theme(axis.text = element_blank(), axis.ticks = element_blank(), plot.margin = margin(0, 3, 0, 
              3))
      lower <- ggplot(data.frame(rank = seq_along(ranks), value = as.numeric(ranks)), aes(rank, value)) + 
          geom_area(fill = "grey70", linewidth = 0) + geom_hline(yintercept = 0, linewidth = 0.25) + scale_x_continuous(limits = c(0, 
          length(ranks)), expand = c(0, 0)) + labs(x = "Rank in ordered dataset", y = "Ranked-list metric") + 
          theme_base
      upper <- upper + theme(plot.title = element_blank())
      header <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = title, size = 8/ggplot2::.pt, fontface = "bold", 
          family = "Arial") + xlim(0, 1) + ylim(0, 1) + theme_void()
      plots[[i]] <- wrap_plots(header, upper, ticks, lower, ncol = 1, heights = c(0.7, 2.8, 0.5, 1.2))
  }

  combined <- wrap_plots(plots, nrow = 1)
  list("Figure S8b | Representative pathway enrichment" = combined)
}

write_s8_svg <- function(plot, path) {
  # Drawing coordinates are in points; physical dimensions determine placement.
  ggplot2::ggsave(path, plot, device = svglite::svglite,
                 width = 487.3 / 72, height = 252.43 / 72, units = "in", bg = "white")
  svg <- paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
  root_match <- regexpr("<svg\\b[^>]*>", svg, perl = TRUE)
  if (root_match[1L] < 0L) stop("SVG root element was not generated.", call. = FALSE)
  root <- regmatches(svg, root_match)
  root <- sub("width=['\"][^'\"]+['\"]", sprintf("width='%.15gmm'", 6191250 / 36000), root)
  root <- sub("height=['\"][^'\"]+['\"]", sprintf("height='%.15gmm'", 3209925 / 36000), root)
  root <- sub(">$", " preserveAspectRatio='none'>", root)
  regmatches(svg, root_match) <- root
  writeLines(enc2utf8(svg), path, useBytes = TRUE)
  invisible(path)
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
  cat("Supplementary Fig. 8b | Representative enrichment curves\nStatistics are recalculated from raw counts and fixed gene sets. Nominal P values reaching the computational floor are reported as <= 1e-10.\n")
  show("representative_pathways", "Representative pathway statistics", c("ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust"))
  cat("\nComplete result tables: ", paste(names(result$tables), collapse = ", "), "\n", sep = "")
  cat("Results directory: ", result$output_dir, "\n", sep = "")
  invisible(result)
}

# Source calculates the results, prints the analysis summary and displays the plots.
if (isTRUE(getOption("SLC3010.autorun", TRUE))) {
  analysis_result <- figure_result <- NULL
  analysis_result <- local({
    needed <- c("digest", "readxl", "DESeq2", "fgsea", "BiocParallel")
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
    paths <- c(frames[basename(frames) == "Fig_S08B_representative_GSEA.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_S08B_representative_GSEA.R"], join(getwd(), "Fig_S08B_representative_GSEA.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_S08B_representative_GSEA.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply one workbook path or directory, or use Source in RStudio.", call. = FALSE)
    path <- path.expand(if (length(args)) args[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR"))
    book <- "Source_Data_RNA_inputs.xlsx"
    sheets <- c("MP_CD8_raw_counts", "S8_gene_ID_mapping", "S8_GSEA_gene_sets", "S8_GSEA_parameters", "S8_pathway_annotations")
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
      run <- join(base, "Fig_S08B_representative_GSEA", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
      needed <- c("ggplot2", "patchwork", "svglite")
      missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
      if (length(missing)) stop("Install these plotting packages: ", paste(missing, collapse = ", "))
      if (interactive() && identical(Sys.getenv("RSTUDIO"), "1") && !identical(names(grDevices::dev.cur()), "RStudioGD"))
        grDevices::dev.new(noRStudioGD = FALSE)
      measure <- !interactive() && grDevices::dev.cur() == 1L
      if (measure) grDevices::pdf(NULL)
      plots <- tryCatch(suppressPackageStartupMessages(plot_results(analysis_result)), finally = if (measure) grDevices::dev.off())
      for (i in seq_along(plots)) {
        size <- analysis_result$sizes[[names(plots)[i]]]
        if (is.null(size)) size <- c(487.3 / 72, 252.43 / 72)
        p <- plots[[i]]
        if (is.function(p)) p <- grid::grid.grabExpr(p(), width = size[1L], height = size[2L])
        if (interactive() && grid::is.grob(p)) { grid::grid.newpage(); grid::grid.draw(p) } else if (interactive()) print(p)
        write_s8_svg(p, file.path(analysis_result$output_dir, sprintf("plot_%02d.svg", i)))
      }
      plots
    }), error = function(e) stop("Analysis tables saved in ", analysis_result$output_dir, ". Plotting failed: ", conditionMessage(e), call. = FALSE))
    figure_result <- analysis_result
  }
}
