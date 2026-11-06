# Supplementary Fig. 4a-b | Bulk-TIL expression programs
# Filter the bulk-TIL count matrix, fit the DESeq2 condition model and calculate gene-wise
# z-scores and program means from rlog-transformed expression.
# Each sample is one mouse. The program scores are descriptive averages of gene-wise z-scores.
#
# Input: put Source_Data_RNA_inputs.xlsx beside this script or in the repository's data/
# folder.
# Numerical packages: digest, readxl / DESeq2, SummarizedExperiment (Bioconductor).
# Plotting packages: dplyr, ggplot2, scales, tidyr.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$individual_program_scores
#   analysis_result$output_dir   # Folder containing Fig_S04B_individual_program_scores.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Supplementary Figure 4a,b: bulk-TIL gene-wise rlog heatmap and complete-program means.
# Inputs: RNA_TIL_counts, RNA_TIL_metadata and RNA_TIL_programs in Source Data; nine independent mice.
# Retains count>=10 in>=3 mice, the original DESeq2/rlog fit, frozen MSigDB 2026.1.Mm memberships,
# highest pooled-rlog-variance gene selection and mean gene-wise z-score calculations.
# Color clipping affects display only; program scores are descriptive means.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  # Excel ranges used for named datasets in this analysis.
  table_ranges <- c(
    "RNA_TIL_counts" = "'TIL_raw_counts'!A1:K45778",
    "RNA_TIL_programs" = "'TIL_program_members'!A1:F9",
    "RNA_TIL_metadata" = "'Sample_information'!A25:O34"
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

  assemble_gene_sets <- function(z) {
    x <- z$table
    order <- unique(x$Term_order)
    sets <- lapply(order, function(i) {
      ix <- which(x$Term_order == i)
      if (length(unique(x$Term[ix])) != 1L || !identical(x$Chunk[ix], seq_along(ix)))
        stop("Inconsistent gene-set chunk order.", call. = FALSE)
      unname(unlist(z$members[ix], use.names = FALSE))
    })
    names(sets) <- vapply(order, function(i) x$Term[match(i, x$Term_order)], character(1))
    if (anyDuplicated(names(sets))) stop("Duplicated gene-set name.", call. = FALSE)
    sets
  }

  read_til_programs <- function() {
    chunks <- read_gene_set_chunks("RNA_TIL_programs")
    sets <- assemble_gene_sets(chunks)
    x <- chunks$table
    if (!"Collection" %in% names(x)) stop("Missing program collection.", call. = FALSE)
    collection <- x$Collection[match(names(sets), x$Term)]
    list(MH = sets[collection == "MH"], GO_targeted = sets[collection == "GO_targeted"])
  }

  suppressPackageStartupMessages(library(DESeq2))
  out <- output_dir

  md_input <- as.data.frame(readxl::read_excel(workbook, sheet = "Sample_information", range = "A25:O34", .name_repair = "minimal"), stringsAsFactors = FALSE)
  md_input$Sample_number <- as.integer(md_input$Sample_number)
  md_input$replicate <- as.integer(md_input$replicate)
  sets <- read_til_programs()
  raw <- as.data.frame(readxl::read_excel(workbook, sheet = "TIL_raw_counts", range = "A1:K45778", .name_repair = "minimal"), stringsAsFactors = FALSE)

  # Map count-matrix columns to the nine mice in the treatment groups.
  input_columns <- c(
    "Beta IL-2 replicate 1 read count" = "Beta IL-2_1",
    "Beta IL-2 replicate 2 read count" = "Beta IL-2_2",
    "Beta IL-2 replicate 3 read count" = "Beta IL-2_3",
    "SLC-3010 replicate 1 read count" = "SLC-3010_1",
    "SLC-3010 replicate 2 read count" = "SLC-3010_2",
    "SLC-3010 replicate 3 read count" = "SLC-3010_3",
    "5344complex replicate 1 read count" = "5344complex_1",
    "5344complex replicate 2 read count" = "5344complex_2",
    "5344complex replicate 3 read count" = "5344complex_3")
  if (!all(names(input_columns) %in% names(raw))) stop("Expected input headers are missing; use the supplied schema-documented input.")
  names(raw)[match(names(input_columns), names(raw))] <- unname(input_columns)

  ann <- raw[,c("Gene_ID","Gene_Symbol")];ann$Gene_ID<-as.character(ann$Gene_ID)
  cnt <- as.matrix(raw[,md_input$Sample])
  stopifnot(nrow(raw)>0L,identical(colnames(cnt),md_input$Sample),all(is.finite(cnt)),all(cnt>=0),all(cnt==round(cnt)),!anyDuplicated(ann$Gene_ID))
  storage.mode(cnt)<-"integer";rownames(cnt)<-ann$Gene_ID
  md<-md_input;rownames(md)<-md$Sample
  md$condition<-factor(md$condition,levels=c("Beta IL-2","SLC-3010","5344complex"))
  keep<-rowSums(cnt>=10)>=3
  if (!any(keep)) stop("No genes meet the prespecified count filter.")
  stopifnot(all(table(md$condition) == 3L))
  # Fit the count model and obtain regularized-log expression for each mouse.
  dds <- DESeq(DESeqDataSetFromMatrix(cnt[keep, ], md, ~condition),
               fitType = "parametric", sfType = "ratio", test = "Wald")
  rld <- rlog(dds, blind = FALSE)
  m<-list(dds=dds,rlog=rld,metadata=md,annotation=ann,sets=sets)

  r <- SummarizedExperiment::assay(m$rlog)
  z <- t(scale(t(r)))
  vars <- apply(r, 1, var)
  IFN <- sort(unique(unlist(sets$MH[c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE")],
      use.names = FALSE)))
  CYCLE <- sort(unique(unlist(sets$MH[c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")], use.names = FALSE)))
  programs <- list(`Immune response\nIFN-alpha/gamma` = IFN, `Proliferation-associated\nE2F/G2M` = CYCLE)
  matrows <- character()
  block <- character()
  for (n in names(programs)) {
      ids <- intersect(programs[[n]], rownames(r))
      ids <- ids[is.finite(vars[ids]) & vars[ids] > 0]
      ids <- ids[order(-vars[ids], ids)]
      matrows <- c(matrows, head(ids, 20))
      block <- c(block, rep(n, min(20, length(ids))))
  }
  wc <- function(x, n) write.csv(x, file.path(out, paste0(n, ".csv")), row.names = FALSE,
      na = "", fileEncoding = "UTF-8")

  gene <- data.frame(program = block, Gene_ID = matrows, Gene_Symbol = ann$Gene_Symbol[match(matrows, ann$Gene_ID)],
      z[matrows, , drop = FALSE], check.names = FALSE)
  wc(gene, "Fig_S04A_gene_z_scores")
  focus <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_E2F_TARGETS",
      "HALLMARK_G2M_CHECKPOINT", "GOBP_T_CELL_ACTIVATION", "GOBP_MITOTIC_CELL_CYCLE")
  lab <- c("IFN-alpha response", "IFN-gamma response", "E2F targets", "G2M checkpoint", "T-cell activation (GO)",
      "Mitotic cell cycle (GO)")
  names(lab) <- focus
  allsets <- c(sets$MH, sets$GO_targeted)
  scores <- do.call(rbind, lapply(focus, function(n) {
      ids <- intersect(allsets[[n]], rownames(z))
      data.frame(md, pathway = n, genes = length(ids), score = colMeans(z[ids, , drop = FALSE]))
  }))
  scores$Program <- unname(lab[scores$pathway])
  wc(scores, "Fig_S04B_individual_program_scores")
  individual <- matrix(NA_real_, length(focus), 9, dimnames = list(focus, md$Sample))
  for (i in seq_len(nrow(scores))) individual[scores$pathway[i], scores$Sample[i]] <- scores$score[i]
  groupmat <- sapply(levels(md$condition), function(k) rowMeans(individual[, md$condition == k, drop = FALSE]))
  gm <- data.frame(pathway = focus, Program = unname(lab), groupmat, check.names = FALSE)
  wc(gm, "Fig_S04B_group_means")
  list(tables = list(group_program_scores = gm, individual_program_scores = scores,
                     gene_z_scores = gene,
                     count_filter = data.frame(input_genes = nrow(cnt), retained_genes = sum(keep), mice = nrow(md))),
       samples = md)
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result) {
  gene <- result$tables$gene_z_scores
  gm <- result$tables$group_program_scores
  md <- result$samples
  # Each column in panel a represents one mouse. The scores in panel b are group means.
  display <- c("Beta IL-2" = "Beta IL-2", "SLC-3010" = "SLC-3010", "5344complex" = "hIL-2/5344 complex")
  gene_data <- tidyr::pivot_longer(gene, cols = dplyr::all_of(md$Sample),
      names_to = "Sample", values_to = "z")
  gene_data$Sample <- factor(gene_data$Sample, levels = md$Sample)
  gene_data$Treatment <- factor(display[as.character(md$condition[match(gene_data$Sample, md$Sample)])],
                                levels = unname(display))
  gene_data$Gene_ID <- factor(gene_data$Gene_ID, levels = rev(unique(gene$Gene_ID)))
  gene_data$program <- factor(gene_data$program, levels = unique(gene$program))
  gene_labels <- setNames(gene$Gene_Symbol, gene$Gene_ID)
  gene_labels[is.na(gene_labels) | !nzchar(gene_labels)] <- names(gene_labels)[is.na(gene_labels) | !nzchar(gene_labels)]
  sample_labels <- setNames(as.character(md$replicate), md$Sample)
  heat_theme <- ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(panel.grid = ggplot2::element_blank(),
      legend.position = "bottom", strip.text.y = ggplot2::element_text(angle = 0),
      plot.caption = ggplot2::element_text(hjust = 0))
  gene_plot <- ggplot2::ggplot(gene_data, ggplot2::aes(Sample, Gene_ID, fill = z)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.16) +
      ggplot2::facet_grid(program ~ Treatment, scales = "free", space = "free") +
      ggplot2::scale_x_discrete(labels = sample_labels) + ggplot2::scale_y_discrete(labels = gene_labels) +
      ggplot2::scale_fill_gradient2("Gene-wise rlog z-score", low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
          midpoint = 0, limits = c(-2, 2), oob = scales::squish, breaks = c(-2, 0, 2),
          labels = c("\u2264 -2", "0", "\u2265 2")) +
      ggplot2::labs(x = "Mouse within treatment", y = NULL,
                    caption = "DESeq2 rlog expression; colour saturates at z = -2 and +2. CSV values retain the full range.") +
      heat_theme
  score_data <- tidyr::pivot_longer(gm, cols = dplyr::all_of(names(display)), names_to = "Treatment", values_to = "Score")
  score_data$Treatment <- factor(display[score_data$Treatment], levels = unname(display))
  score_data$Program <- factor(score_data$Program, levels = rev(gm$Program))
  score_plot <- ggplot2::ggplot(score_data, ggplot2::aes(Treatment, Program, fill = Score)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", Score)), colour = "black", size = 3.3) +
      ggplot2::scale_fill_gradient2("Mean gene-wise rlog z-score", low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
          midpoint = 0, limits = c(-1, 1), oob = scales::squish, breaks = c(-1, 0, 1), labels = c("\u2264 -1", "0", "\u2265 1")) +
      ggplot2::labs(x = NULL, y = NULL, caption = "Each tile summarizes the three mice in its treatment group.") +
      heat_theme
  list("Figure S4a | Individual gene expression" = gene_plot,
       "Figure S4b | Tumor gene programs" = score_plot)
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
  cat("Supplementary Fig. 4a-b | Bulk-TIL expression programs\nProgram scores are descriptive averages of gene-wise z-scores from rlog-transformed expression; each sample is one mouse.\n")
  show("count_filter", "Count filtering and biological sample size", c("input_genes", "retained_genes", "mice"))
  show("group_program_scores", "Group means of expression-program scores", c("Program", "Beta IL-2", "SLC-3010", "5344complex"))
  cat("\nComplete result tables: ", paste(names(result$tables), collapse = ", "), "\n", sep = "")
  cat("Results directory: ", result$output_dir, "\n", sep = "")
  invisible(result)
}

# Source calculates the results, prints the analysis summary and displays the plots.
if (isTRUE(getOption("SLC3010.autorun", TRUE))) {
  analysis_result <- figure_result <- NULL
  analysis_result <- local({
    needed <- c("DESeq2", "SummarizedExperiment", "digest", "readxl")
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
    paths <- c(frames[basename(frames) == "Fig_S04A-B_bulk_TIL_transcriptome.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_S04A-B_bulk_TIL_transcriptome.R"], join(getwd(), "Fig_S04A-B_bulk_TIL_transcriptome.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_S04A-B_bulk_TIL_transcriptome.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply one workbook path or directory, or use Source in RStudio.", call. = FALSE)
    path <- path.expand(if (length(args)) args[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR"))
    book <- "Source_Data_RNA_inputs.xlsx"
    sheets <- c("Sample_information", "TIL_program_members", "TIL_raw_counts")
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
      run <- join(base, "Fig_S04A-B_bulk_TIL_transcriptome", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
      needed <- c("dplyr", "ggplot2", "scales", "tidyr")
      missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
      if (length(missing)) stop("Install these plotting packages: ", paste(missing, collapse = ", "))
      if (interactive() && identical(Sys.getenv("RSTUDIO"), "1") && !identical(names(grDevices::dev.cur()), "RStudioGD"))
        grDevices::dev.new(noRStudioGD = FALSE)
      measure <- !interactive() && grDevices::dev.cur() == 1L
      if (measure) grDevices::pdf(NULL)
      plots <- tryCatch(suppressPackageStartupMessages(plot_results(analysis_result)), finally = if (measure) grDevices::dev.off())
      for (i in seq_along(plots)) {
        size <- analysis_result$sizes[[names(plots)[i]]]
        if (is.null(size)) size <- c(10, 8)
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
