# Fig. 5h | MP-CD8 principal component analysis
# Calculate centered, unscaled PCA on log2-transformed normalized counts for up to 5,000
# variable genes.
# PCA uses the supplied normalized matrix. A separate check recalculates normalization from
# raw counts across all 18 samples.
#
# Input: put Source_Data_RNA_inputs.xlsx beside this script or in the repository's data/
# folder.
# Numerical packages: digest, dplyr, readr, readxl, stringr, tibble / DESeq2 (Bioconductor).
# Plotting packages: ggplot2, ggrepel, purrr, tibble.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$pca_coordinates
#   analysis_result$output_dir   # Folder containing Fig_05H_PCA_coordinates.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Figure 5h: PCA of the 12 RNA-seq samples.
# Inputs: RNA_PCA_counts, RNA_MP_raw_counts and RNA_MP_samples in the RNA Source Data workbook.
# PCA uses the supplied normalized counts for 12 samples; an RLE calculation checks normalization across all 18 samples.
# Duplicate symbols are averaged; missing values become zero; log2(count+1), up to 5,000 variable genes, centered unscaled PCA.
# Shows sample labels and descriptive ellipses by sampling time and genotype.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  # Excel ranges used for named datasets in this analysis.
  table_ranges <- c(
    "RNA_PCA_counts" = "'MP_CD8_normalized_counts'!A1:T13137",
    "RNA_MP_raw_counts" = "'MP_CD8_raw_counts'!A1:T45775",
    "RNA_MP_samples" = "'Sample_information'!A3:D21"
  )
  INPUT_DATASET <- "RNA_PCA_counts"

  USE_LOG2_TRANSFORM <- TRUE

  PSEUDOCOUNT <- 1

  N_TOP_VAR_GENES <- 5000

  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(stringr))
  OUT_DIR <- output_dir

  require_columns <- function(x, cols, filename) {
      missing <- setdiff(cols, names(x))
      if (length(missing)) {
          stop("Missing columns in ", filename, ": ", paste(missing, collapse = ", "), "\nExisting columns: ", paste(names(x), 
              collapse = ", "), call. = FALSE)
      }
  }

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

  # Read and normalize the MP-CD8 expression measurements.
  expr_data <- readxl::read_excel(workbook, range = table_ranges[[INPUT_DATASET]], .name_repair = "minimal")

  # Recalculate RLE normalization across all 18 samples and compare it with the supplied matrix.
  compare_normalization <- function(workbook, expr_data, output_dir) {
    # Keep the caller's RNG state, including the state with no seed yet.
    had_rng_state <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    previous_rng_state <- if (had_rng_state) get(".Random.seed", envir = .GlobalEnv) else NULL
    previous_rng_kind <- RNGkind()
    on.exit({
      do.call(RNGkind, as.list(previous_rng_kind))
      if (had_rng_state) {
        assign(".Random.seed", previous_rng_state, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(list = ".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    write_method_record <- function(output_dir, prefix, notes) {
      writeLines(c(notes, "", capture.output(utils::sessionInfo())),
                 file.path(output_dir, paste0(prefix, "_method_and_runtime.txt")), useBytes = TRUE)
    }

    comparison_output_dir <- file.path(output_dir, "reproducibility")
    dir.create(comparison_output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(comparison_output_dir))
      stop("Cannot create the reproducibility output directory.", call. = FALSE)

    raw <- as.data.frame(readxl::read_excel(workbook, sheet = "MP_CD8_raw_counts", range = "A1:T45775", .name_repair = "minimal"), stringsAsFactors = FALSE)
    # Sample metadata: Count_matrix_column selects raw-count columns; Sample.ID links
    # normalized counts. Genotype and Time_point describe the corresponding mouse sample.
    samples <- as.data.frame(readxl::read_excel(workbook, sheet = "Sample_information", range = "A3:D21", .name_repair = "minimal"), stringsAsFactors = FALSE)
    count_columns <- as.character(samples$Count_matrix_column)
    if (!all(count_columns %in% names(raw))) stop("Sample key does not match raw count columns.")
    counts <- as.matrix(raw[, count_columns])
    if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0) || any(counts != round(counts))) stop("Invalid raw counts.")
    storage.mode(counts) <- "integer"
    keep <- rowSums(counts == 0L) == 0L
    retained <- counts[keep, , drop = FALSE]
    rownames(retained) <- as.character(raw$Gene_ID[keep])
    size_factors <- DESeq2::estimateSizeFactorsForMatrix(retained)
    normalized <- sweep(retained, 2, size_factors, "/")
    reference_counts <- as.data.frame(expr_data, stringsAsFactors = FALSE)
    columns <- c("wild-type Baseline sample 1 RLE-normalized count", "wild-type Baseline sample 2 RLE-normalized count",
      "wild-type Baseline sample 3 RLE-normalized count", "Triple hIL-2 Baseline sample 4 RLE-normalized count",
      "Triple hIL-2 Baseline sample 5 RLE-normalized count", "Triple hIL-2 Baseline sample 6 RLE-normalized count",
      "wild-type C1D4 sample 1 RLE-normalized count", "wild-type C1D4 sample 2 RLE-normalized count",
      "wild-type C1D4 sample 3 RLE-normalized count", "Triple hIL-2 C1D4 sample 4 RLE-normalized count",
      "Triple hIL-2 C1D4 sample 5 RLE-normalized count", "Triple hIL-2 C1D4 sample 6 RLE-normalized count",
      "wild-type C3D4 sample 1 RLE-normalized count", "wild-type C3D4 sample 2 RLE-normalized count",
      "wild-type C3D4 sample 3 RLE-normalized count", "Triple hIL-2 C3D4 sample 4 RLE-normalized count",
      "Triple hIL-2 C3D4 sample 5 RLE-normalized count", "Triple hIL-2 C3D4 sample 6 RLE-normalized count")
    sample_order <- c(paste0("Baseline_", 1:6), paste0("C1D4_", 1:6), paste0("C3D4_", 1:6))
    names(columns) <- sample_order
    if (!setequal(samples$Sample.ID, sample_order)) stop("Unexpected sample identifiers.")
    reference_matrix <- as.matrix(reference_counts[, unname(columns[as.character(samples$Sample.ID)])])
    matched <- match(rownames(normalized), as.character(reference_counts$Gene_ID))
    shared <- !is.na(matched)
    if (!any(shared)) stop("No common genes for the normalization comparison.")
    if (any(!shared) || sum(shared) != nrow(reference_counts))
      warning("The supplied and recomputed matrices have different genes; comparing ", sum(shared), " common genes.")
    reference_matrix <- reference_matrix[matched[shared], , drop = FALSE]
    difference <- normalized[shared, , drop = FALSE] - reference_matrix
    summary <- data.frame(Sample = samples$Sample.ID, raw_count_column = count_columns,
      retained_genes = nrow(normalized), size_factor = as.numeric(size_factors),
      exact_matches = colSums(difference == 0), max_absolute_difference = apply(abs(difference), 2, max),
      max_relative_difference = apply(abs(difference)/pmax(abs(reference_matrix), .Machine$double.xmin), 2, max))
    write.csv(summary, file.path(comparison_output_dir, "MP_CD8_normalization_comparison.csv"), row.names = FALSE)
    output <- data.frame(Gene_ID = rownames(normalized), normalized, check.names = FALSE)
    names(output)[-1] <- as.character(samples$Sample.ID)
    write.csv(output, file.path(comparison_output_dir, "MP_CD8_reestimated_normalized_counts.csv"), row.names = FALSE)
    write_method_record(comparison_output_dir, "MP_CD8", c(
      "Independent comparison of the documented preprocessing with the supplied RLE-normalized count table.",
      "Gene filter: remove any gene with zero reads in any of the 18 samples. No sample filtering.",
      "Normalization: current DESeq2 estimateSizeFactorsForMatrix, using its standard median-ratio/RLE method.",
      "The original DESeq2 version, fitting design and contrast-fitting scope are not established by the supplied report.",
      "Figure 5h uses the supplied normalized counts; this recalculation is a normalization comparison."))
    summary
  }
  normalization <- compare_normalization(workbook, expr_data, output_dir)
  # Read the normalized expression matrix for the PCA.
  input_columns <- c(
    "wild-type Baseline sample 1 RLE-normalized count" = "N_Baseline_1",
    "wild-type Baseline sample 2 RLE-normalized count" = "N_Baseline_2",
    "wild-type Baseline sample 3 RLE-normalized count" = "N_Baseline_3",
    "wild-type C1D4 sample 1 RLE-normalized count" = "N_C1D4_1",
    "wild-type C1D4 sample 2 RLE-normalized count" = "N_C1D4_2",
    "wild-type C1D4 sample 3 RLE-normalized count" = "N_C1D4_3",
    "wild-type C3D4 sample 1 RLE-normalized count" = "N_C3D4_1",
    "wild-type C3D4 sample 2 RLE-normalized count" = "N_C3D4_2",
    "wild-type C3D4 sample 3 RLE-normalized count" = "N_C3D4_3",
    "Triple hIL-2 Baseline sample 4 RLE-normalized count" = "N_Baseline_4",
    "Triple hIL-2 Baseline sample 5 RLE-normalized count" = "N_Baseline_5",
    "Triple hIL-2 Baseline sample 6 RLE-normalized count" = "N_Baseline_6",
    "Triple hIL-2 C1D4 sample 4 RLE-normalized count" = "N_C1D4_4",
    "Triple hIL-2 C1D4 sample 5 RLE-normalized count" = "N_C1D4_5",
    "Triple hIL-2 C1D4 sample 6 RLE-normalized count" = "N_C1D4_6",
    "Triple hIL-2 C3D4 sample 4 RLE-normalized count" = "N_C3D4_4",
    "Triple hIL-2 C3D4 sample 5 RLE-normalized count" = "N_C3D4_5",
    "Triple hIL-2 C3D4 sample 6 RLE-normalized count" = "N_C3D4_6")
  if (!all(names(input_columns) %in% names(expr_data))) stop("Expected input headers are missing; use the supplied schema-documented input.")
  names(expr_data)[match(names(input_columns), names(expr_data))] <- unname(input_columns)


  require_columns(expr_data, c("Gene_ID", "Gene_Symbol"), INPUT_DATASET)

  sample_cols <- names(expr_data)[stringr::str_detect(names(expr_data), "Baseline|C1D4")]

  if (length(sample_cols) < 3L) stop("Fewer than three Baseline/C1D4 sample columns.")

  safe_gene_label <- function(gene_symbol, gene_id) {
      gene_symbol <- as.character(gene_symbol)
      gene_id <- as.character(gene_id)
      ifelse(is.na(gene_symbol) | gene_symbol == "" | gene_symbol == "NA", gene_id, gene_symbol)
  }

  expr_tbl <- expr_data %>% mutate(Gene_Label = safe_gene_label(Gene_Symbol, Gene_ID)) %>% select(Gene_Label, all_of(sample_cols)) %>% 
      mutate(across(all_of(sample_cols), ~suppressWarnings(as.numeric(.x)))) %>% group_by(Gene_Label) %>% summarise(across(all_of(sample_cols), 
      ~mean(.x, na.rm = TRUE)), .groups = "drop")

  if (any(is.na(expr_tbl$Gene_Label) | expr_tbl$Gene_Label == "")) {
      stop("Some genes have neither a usable Gene_Symbol nor Gene_ID.")
  }

  expr_mat <- as.matrix(select(expr_tbl, all_of(sample_cols)))

  rownames(expr_mat) <- expr_tbl$Gene_Label

  n_nonfinite <- sum(!is.finite(expr_mat))

  expr_mat[!is.finite(expr_mat)] <- NA

  expr_mat[is.na(expr_mat)] <- 0

  if (n_nonfinite) warning(n_nonfinite, " non-finite expression values replaced by 0.")

  if (any(expr_mat < 0)) stop("Negative normalized counts found; inspect the input before log2 transform.")

  if (USE_LOG2_TRANSFORM) expr_mat <- log2(expr_mat + PSEUDOCOUNT)

  vars <- apply(expr_mat, 1, stats::var, na.rm = TRUE)

  expr_mat <- expr_mat[is.finite(vars) & vars > 0, , drop = FALSE]

  if (is.finite(N_TOP_VAR_GENES) && nrow(expr_mat) > N_TOP_VAR_GENES) {
      vars <- apply(expr_mat, 1, stats::var, na.rm = TRUE)
      expr_mat <- expr_mat[order(vars, decreasing = TRUE)[seq_len(N_TOP_VAR_GENES)], , drop = FALSE]
  }

  if (nrow(expr_mat) < 2L) stop("Not enough non-zero-variance genes for PCA.")

  expr_t <- t(expr_mat)

  # Link displayed samples to the biological identifiers and metadata supplied with the counts.
  metadata <- as.data.frame(readxl::read_excel(workbook, range = table_ranges[["RNA_MP_samples"]],
      .name_repair = "minimal"), stringsAsFactors = FALSE)
  require_columns(metadata, c("Sample.ID", "Genotype", "Time_point"), "MP-CD8 sample information")
  sample_ids <- stringr::str_remove(rownames(expr_t), "^N[_-]?")
  if (anyDuplicated(metadata$Sample.ID) || anyDuplicated(sample_ids)) stop("Duplicated biological sample identifiers.")
  positions <- match(sample_ids, metadata$Sample.ID)
  if (anyNA(positions)) stop("PCA samples are missing from Sample_information.")
  sample_info <- tibble::tibble(Sample.ID = sample_ids,
      Genotype = as.character(metadata$Genotype[positions]),
      Time_point = as.character(metadata$Time_point[positions]))
  if (anyNA(sample_info) || !all(sample_info$Genotype %in% c("Wild type", "Triple hIL-2")) ||
      !all(sample_info$Time_point %in% c("Baseline", "C1D4")))
      stop("Unrecognized PCA genotype or sampling time in Sample_information.")
  # Plotting keys retain the established palette, ordering and ellipse definitions.
  sample_info <- sample_info %>% mutate(Group = Time_point,
      Strain = if_else(Genotype == "Wild type", "wt", "Triple hIL-2"),
      Group2 = paste(Group, Strain, sep = "_"),
      Group3 = case_when(Group == "Baseline" ~ "Baseline",
          Group == "C1D4" & Strain == "wt" ~ "C1D4_wt",
          Group == "C1D4" & Strain == "Triple hIL-2" ~ "C1D4_Triple hIL-2"),
      Label = Sample.ID)

  group_levels <- c("Baseline_wt", "Baseline_Triple hIL-2", "C1D4_wt", "C1D4_Triple hIL-2")

  sample_info$Group2 <- factor(sample_info$Group2, levels = group_levels)

  sample_info$Group3 <- factor(sample_info$Group3, levels = c("Baseline", "C1D4_wt", "C1D4_Triple hIL-2"))

  pca_res <- stats::prcomp(expr_t, center = TRUE, scale. = FALSE)

  pca_var <- round(100 * pca_res$sdev^2/sum(pca_res$sdev^2), 1)

  pca_df <- bind_cols(tibble::as_tibble(pca_res$x[, 1:2, drop = FALSE]), sample_info)

  variance <- data.frame(PC = seq_along(pca_var), variance_percent = pca_var)
  coordinates <- as.data.frame(select(pca_df, Sample.ID, Genotype, Time_point, PC1, PC2))
  write_table(coordinates, "Fig_05H_PCA_coordinates.csv")
  write_table(variance, "Fig_05H_PCA_variance.csv")
  list(tables = list(pca_variance = variance, pca_coordinates = coordinates,
                     normalization_comparison = normalization),
       plot_data = pca_df, model = pca_res)
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result) {
  suppressPackageStartupMessages(library(ggplot2))
  pca_df <- result$plot_data
  pca_var <- result$tables$pca_variance$variance_percent
  SHOW_ELLIPSES <- TRUE

  ELLIPSE_LEVELS <- c(Baseline = 0.8, C1D4_wt = 0.8, `C1D4_Triple hIL-2` = 0.68)
  make_cov_ellipse <- function(df, group_name, level = 0.8, n = 120) {
      if (nrow(df) < 3L || stats::sd(df$PC1) == 0 || stats::sd(df$PC2) == 0) 
          return(NULL)
      cov_mat <- stats::cov(cbind(df$PC1, df$PC2), use = "complete.obs")
      if (any(!is.finite(cov_mat))) 
          return(NULL)
      eig <- tryCatch(eigen(cov_mat), error = function(e) NULL)
      if (is.null(eig) || any(eig$values <= 0)) 
          return(NULL)
      theta <- seq(0, 2 * pi, length.out = n)
      circle <- cbind(cos(theta), sin(theta))
      ell <- sweep(circle %*% diag(sqrt(eig$values)) %*% t(eig$vectors) * sqrt(stats::qchisq(level, df = 2)), 2, 
          c(mean(df$PC1), mean(df$PC2)), "+")
      tibble::tibble(PC1 = ell[, 1], PC2 = ell[, 2], EllipseGroup = group_name)
  }

  ellipse_df <- tibble::tibble()

  if (SHOW_ELLIPSES) {
      ellipse_df <- purrr::map_dfr(split(pca_df, pca_df$Group3), function(d) {
          if (!nrow(d)) 
              return(NULL)
          g <- as.character(d$Group3[1])
          make_cov_ellipse(d, g, ELLIPSE_LEVELS[[g]])
      })
  }

  group2_colors <- c(Baseline_wt = "#D0D0D0", `Baseline_Triple hIL-2` = "#222222", C1D4_wt = "#E69F00", `C1D4_Triple hIL-2` = "#E41A1C")

  ellipse_colors <- c(Baseline = "#555555", C1D4_wt = "#E69F00", `C1D4_Triple hIL-2` = "#E41A1C")

  group2_labels <- c("Baseline wt", "Baseline triple hIL-2", "C1D4 wt", "C1D4 triple hIL-2")

  xx <- c(pca_df$PC1, if (nrow(ellipse_df)) ellipse_df$PC1 else numeric())

  yy <- c(pca_df$PC2, if (nrow(ellipse_df)) ellipse_df$PC2 else numeric())

  pad_range <- function(z, fraction = 0.12) {
      r <- range(z, finite = TRUE)
      pad <- max(diff(r) * fraction, 0.5)
      r + c(-pad, pad)
  }

  xlim_pca <- pad_range(xx)

  ylim_pca <- pad_range(yy)

  pca_plot <- ggplot(pca_df, aes(PC1, PC2))
  if (SHOW_ELLIPSES && nrow(ellipse_df)) {
      pca_plot <- pca_plot + geom_path(data = ellipse_df, aes(PC1, PC2, colour = EllipseGroup),
          inherit.aes = FALSE, linewidth = 0.45, linetype = "22", show.legend = FALSE)
  }
  pca_plot <- pca_plot +
      geom_point(aes(fill = Group2), colour = "black", shape = 21, size = 3, stroke = 0.35) +
      ggrepel::geom_text_repel(aes(label = Label), size = 3, seed = 123, max.overlaps = Inf,
          min.segment.length = 0, segment.colour = "grey55", max.iter = 3000) +
      scale_colour_manual(values = ellipse_colors, guide = "none") +
      scale_fill_manual(name = NULL, values = group2_colors, breaks = names(group2_colors),
                        labels = group2_labels, drop = FALSE) +
      coord_cartesian(xlim = xlim_pca, ylim = ylim_pca) +
      labs(x = sprintf("PC1 (%s%%)", pca_var[1]), y = sprintf("PC2 (%s%%)", pca_var[2]),
           caption = "Descriptive ellipses: pooled Baseline and C1D4 WT, 80%; C1D4 triple hIL-2, 68%.") +
      theme_classic(base_size = 10) + theme(legend.position = "bottom", plot.caption = element_text(hjust = 0)) +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE))
  list("Figure 5h | MP-CD8 PCA" = pca_plot)
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
  cat("Fig. 5h | MP-CD8 principal component analysis\nCentered, unscaled PCA uses log2-normalized counts. Normalization covers all 18 samples; this display includes the 12 Baseline/C1D4 samples.\n")
  show("pca_variance", "Variance explained by each principal component", c("PC", "variance_percent"))
  show("pca_coordinates", "Coordinates of the 12 displayed samples", c("Sample.ID", "Genotype", "Time_point", "PC1", "PC2"))
  cat("\nComplete result tables: ", paste(names(result$tables), collapse = ", "), "\n", sep = "")
  cat("Results directory: ", result$output_dir, "\n", sep = "")
  invisible(result)
}

# Source calculates the results, prints the analysis summary and displays the plots.
if (isTRUE(getOption("SLC3010.autorun", TRUE))) {
  analysis_result <- figure_result <- NULL
  analysis_result <- local({
    needed <- c("DESeq2", "digest", "dplyr", "readr", "readxl", "stringr", "tibble")
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
    paths <- c(frames[basename(frames) == "Fig_05H_PCA.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_05H_PCA.R"], join(getwd(), "Fig_05H_PCA.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_05H_PCA.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply one workbook path or directory, or use Source in RStudio.", call. = FALSE)
    path <- path.expand(if (length(args)) args[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR"))
    book <- "Source_Data_RNA_inputs.xlsx"
    sheets <- c("MP_CD8_normalized_counts", "MP_CD8_raw_counts", "Sample_information")
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
      run <- join(base, "Fig_05H_PCA", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
      needed <- c("ggplot2", "ggrepel", "purrr", "tibble")
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
