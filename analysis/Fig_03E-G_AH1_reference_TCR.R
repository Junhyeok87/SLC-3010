# Fig. 3e-g | AH1-reference TCR sequences detected per tumor
# Sum eligible TRB sequence counts by CDR3 sequence within each tumor; compare observed match
# counts between groups with an unpaired two-sided Student's t-test with equal variances.
# The Student's t-test uses observed richness. Expected richness at a common sequencing depth is a
# separate per-tumor estimate.
#
# Input: put Source_Data_TCR_inputs.xlsx beside this script or in the repository's data/
# folder.
# Numerical packages: digest, readxl.
# Plotting packages: ggplot2, ggrepel, scales.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$student_comparison
#   analysis_result$output_dir   # Folder containing Fig_03E_Student_contrast.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Figure 3e–g: AH1-reference TCR detection and sequencing-depth comparison.
# Clone observations, sample metadata and reference CDR3 sequences are read from Source Data Excel.
# Counts are aggregated by eligible TRB CDR3 sequence. Common-depth richness uses the
# hypergeometric detection expectation; panel e uses an unadjusted unpaired two-sided Student's t-test with equal variances.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  # Clone-table bounds are keyed by the worksheet names recorded for each sample.
  clone_ranges <- c(
    "V1_clone_counts" = "A1:G24921",
    "V2_clone_counts" = "A1:G20360",
    "V3_clone_counts" = "A1:G13422",
    "V4_clone_counts" = "A1:G12519",
    "S1_clone_counts" = "A1:G83152",
    "S2_clone_counts" = "A1:G52598",
    "S3_clone_counts" = "A1:G61038"
  )
  # Calculate the sequence-matching and per-tumor statistics.
  summarize_tcr <- function(summary, M, totals, expected, per_sequence, out) {
    frac <- sweep(M,2,totals,"/")
    thresholded <- M>=2 & frac>=1e-5
    summary$AH1_Exact_CDR3_Matches_Any <- colSums(M>0)
    v <- summary$AH1_Exact_CDR3_Matches_Any[summary$Group=="Vehicle"]
    z <- summary$AH1_Exact_CDR3_Matches_Any[summary$Group=="SLC-3010"]
    test <- stats::t.test(v, z, alternative = "two.sided", paired = FALSE, var.equal = TRUE)
    pv <- test$p.value
    statistics <- data.frame(Endpoint="AH1_Exact_CDR3_Matches_Any",Vehicle_n=length(v),SLC3010_n=length(z),Vehicle_mean=mean(v),SLC3010_mean=mean(z),P_Student_two_sided=pv)
    grid <- do.call(rbind,per_sequence)
    grid$Sample_label <- expected$Sample_label[match(grid$Sample,expected$Sample)]
    grid$TRB_ppm <- 1e6*grid$TRB_fraction
    grid$Heat_TRB <- log10(grid$TRB_ppm+1)
    detected <- rowSums(M>0)>0
    ord <- order(-rowSums(thresholded),-rowSums(frac),rownames(M),method="radix")
    order_all <- rownames(M)[ord]
    order_detected <- order_all[order_all %in% rownames(M)[detected]]
    utils::write.csv(summary,file.path(out,"Fig_03E-F_per_tumor_results.csv"),row.names=FALSE,na="")
    utils::write.csv(statistics,file.path(out,"Fig_03E_Student_statistics.csv"),row.names=FALSE,na="")
    utils::write.csv(grid,file.path(out,"Fig_03G_sequence_counts_frequencies.csv"),row.names=FALSE,na="")
  contrast <- data.frame(
        group1 = "Vehicle", group2 = "SLC-3010", n1 = length(v), n2 = length(z),
        mean1 = mean(v), mean2 = mean(z), vehicle_minus_slc = mean(v) - mean(z),
        ci95_low = test$conf.int[1], ci95_high = test$conf.int[2],
        t = unname(test$statistic), df = unname(test$parameter), p = pv)
      utils::write.csv(contrast, file.path(out, "Fig_03E_Student_contrast.csv"), row.names = FALSE)
      list(summary = summary, contrast = contrast, frequencies = grid, detected_cdr3 = order_detected)
  }
  COUNT_COLUMN <- "Sequence_Count"
  MIN_CDR3_AA <- 8L
  MAX_CDR3_AA <- 20L
  fail <- function(...) stop(paste0(...), call. = FALSE)
  assert <- function(ok, ...) if (!isTRUE(ok)) fail(...)
  clean <- function(x) trimws(as.character(x))
  aa_clean <- function(x) {
    value <- toupper(clean(x));value[is.na(x) | !nzchar(value)] <- NA_character_;value
  }
  numbers <- function(value,name,integer=FALSE,missing_ok=FALSE) {
    out <- suppressWarnings(as.numeric(value))
    missing <- is.na(value) | !nzchar(trimws(as.character(value)))
    bad <- !missing & (!is.finite(out) | out < 0)
    if(integer)bad <- bad | (!missing & out!=floor(out))
    if(!missing_ok)bad <- bad | missing
    assert(!any(bad),"Invalid ",name," values.")
    out[missing] <- NA_real_;out
  }

  # The deposited clone tables have no explicit productivity annotation.
  # Sequence_Count is the number of collapsed sequences per clonal group.
  # Count-weighted frequencies use its sum across eligible TRB rows in each sample.
  canonicalize <- function(z,sample) {
    required <- c("Clone_ID","Sequence_Count","Copy_Count","Chain","CDR3_aa","V_Gene","J_Gene")
    assert(all(required %in% names(z)) && nrow(z)>0L,"Invalid clone table for ",sample)
    chain <- toupper(clean(z$Chain));aa <- aa_clean(z$CDR3_aa)
    count <- numbers(z[[COUNT_COLUMN]],paste0(sample," ",COUNT_COLUMN),integer=TRUE,missing_ok=TRUE)
    trb <- !is.na(chain) & chain=="TRB"
    standard <- !is.na(aa) & grepl("^[ACDEFGHIKLMNPQRSTVWY]+$",aa)
    inlength <- !is.na(aa) & nchar(aa)>=MIN_CDR3_AA & nchar(aa)<=MAX_CDR3_AA
    positive <- is.finite(count) & count>0
    eligible <- trb & standard & inlength & positive
    assert(any(eligible),"No eligible TRB rows for ",sample)
    value <- count[eligible];sequence <- aa[eligible]
    assert(sum(value)<2^53,"Count total exceeds exact integer precision.")
    combined <- rowsum(value,group=sequence,reorder=TRUE)
    tab <- data.frame(CDR3_aa=rownames(combined),Count=as.numeric(combined[,1L]),stringsAsFactors=FALSE)
    tab <- tab[order(tab$CDR3_aa,method="radix"),,drop=FALSE];rownames(tab)<-NULL
    assert(sum(tab$Count)==sum(value),"Count total changed during CDR3 aggregation.")
    qc <- data.frame(Sample=sample,Raw_rows=nrow(z),TRB_rows=sum(trb),
      Nonproductive_TRB_rows=NA_integer_,Invalid_AA_TRB_rows=sum(trb & !standard),
      Outside_length_TRB_rows=sum(trb & standard & !inlength),
      Missing_or_zero_count_TRB_rows=sum(trb & !positive),Eligible_TRB_rows=sum(eligible),
      Distinct_TRB_sequences=nrow(tab),TRB_Count=sum(value),Count_basis=COUNT_COLUMN,
      CDR3_min_aa=MIN_CDR3_AA,CDR3_max_aa=MAX_CDR3_AA,Productivity_fields="",
      Productivity_assessed=FALSE,
      Eligible_rows_without_explicit_productivity=sum(eligible),stringsAsFactors=FALSE)
    list(counts=tab,qc=qc)
  }
  expected_breadth <- function(matched_counts,N,depth) {
    assert(length(N)==1L && is.finite(N) && N>0 && N==floor(N) && N<2^53,"Invalid total TRB depth.")
    assert(length(depth)==1L && is.finite(depth) && depth>=0 && depth<=N && depth==floor(depth),"Invalid common depth.")
    assert(all(is.finite(matched_counts) & matched_counts>=0 & matched_counts==floor(matched_counts)) && sum(matched_counts)<=N,"Invalid reference counts.")
    positive <- matched_counts[matched_counts>0]
    if(!length(positive) || depth==0)return(0)
    sum(stats::phyper(0,m=positive,n=N-positive,k=depth,lower.tail=FALSE))
  }

  samples <- as.data.frame(readxl::read_excel(workbook, sheet = "Sample_information", range = "A3:F10", col_types="text", .name_repair = "minimal"))
  assert(identical(names(samples), c("Sample", "Group", "Replicate", "Sample_label",
                                    "CD8_Tcells_million_per_g", "Clonotype_table")),
         "Unexpected TCR sample metadata columns.")
  assert(nrow(samples)==7L && !anyDuplicated(samples$Sample),"Seven distinct tumor samples are required.")
  assert(!anyNA(samples$Clonotype_table) && !anyDuplicated(samples$Clonotype_table) &&
         setequal(samples$Clonotype_table, names(clone_ranges)),
         "Seven distinct, recognized clonotype worksheets are required.")
  expected_samples <- samples$Sample
  expected_groups <- samples$Group
  expected_labels <- samples$Sample_label
  assert(sum(expected_groups=="Vehicle")==4L && sum(expected_groups=="SLC-3010")==3L,"Unexpected treatment groups.")
  expected <- samples[,c("Sample","Group","Replicate","Sample_label"),drop=FALSE]
  flow <- samples[,c("Sample","Group","CD8_Tcells_million_per_g"),drop=FALSE]
  flow$CD8_Tcells_million_per_g <- numbers(flow$CD8_Tcells_million_per_g,"CD8 T-cell density")
  canonical <- lapply(seq_len(nrow(samples)),function(i) {
    clone_sheet <- samples$Clonotype_table[i]
    table <- as.data.frame(readxl::read_excel(workbook, sheet = clone_sheet,
                           range = clone_ranges[[clone_sheet]], col_types="text", .name_repair = "minimal"))
    canonicalize(table,samples$Sample[i])
  })
  totals <- vapply(canonical,function(x)sum(x$counts$Count),numeric(1))
  common_depth <- min(totals)
  ref <- as.data.frame(readxl::read_excel(workbook, sheet = "AH1_reference", range = "A1:B41", col_types="text", .name_repair = "minimal"))
  ref$CDR3_aa <- aa_clean(ref$CDR3_aa)
  assert(nrow(ref)==40L && !anyDuplicated(ref$CDR3_aa) && all(!is.na(ref$CDR3_aa) & grepl("^[ACDEFGHIKLMNPQRSTVWY]+$",ref$CDR3_aa)),"The AH1 reference must contain 40 distinct CDR3 sequences.")
  ref <- ref[order(ref$CDR3_aa,method="radix"),,drop=FALSE];rownames(ref)<-NULL
  M <- matrix(0,nrow=nrow(ref),ncol=7L,dimnames=list(ref$CDR3_aa,expected_samples))
  per_sequence <- vector("list",7L)
  for(i in seq_len(7L)) {
    tab <- canonical[[i]]$counts;position <- match(ref$CDR3_aa,tab$CDR3_aa)
    hit <- !is.na(position);M[hit,i] <- tab$Count[position[hit]]
    probability <- stats::phyper(0,m=M[,i],n=totals[i]-M[,i],k=common_depth,lower.tail=FALSE)
    per_sequence[[i]] <- data.frame(Sample=expected_samples[i],Group=expected_groups[i],CDR3_aa=ref$CDR3_aa,
      Count=M[,i],TRB_Count=totals[i],Common_TRB_depth=common_depth,
      TRB_fraction=M[,i]/totals[i],Detection_probability=probability,stringsAsFactors=FALSE)
  }
  y <- vapply(seq_len(7L),function(i)expected_breadth(M[,i],totals[i],common_depth),numeric(1))
  observed <- colSums(M>0)
  assert(all(is.finite(y) & y>=0 & y<=observed+1e-10),"Invalid expected richness.")
  assert(all(abs(y[totals==common_depth]-observed[totals==common_depth])<1e-10),"Minimum-depth richness is inconsistent.")
  summary <- data.frame(expected,CD8_Tcells_million_per_g=flow$CD8_Tcells_million_per_g,
    TRB_Count=totals,Common_TRB_depth=common_depth,AH1_Count=colSums(M),AH1_richness=observed,
    AH1_expected_richness_at_common_depth=y,Count_basis=COUNT_COLUMN,stringsAsFactors=FALSE)
  counted <- summarize_tcr(summary, M, totals, expected, per_sequence, output_dir)
  qc <- do.call(rbind, lapply(canonical, function(sample) sample$qc))
  utils::write.csv(qc, file.path(output_dir, "Fig_03E-G_TRB_count_QC.csv"), row.names = FALSE)
  list(tables = list(student_comparison = counted$contrast, per_tumor = counted$summary,
                     count_qc = qc, sequence_frequencies = counted$frequencies),
       detected_cdr3 = counted$detected_cdr3)
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result) {
  count <- result$tables$per_tumor
  common_depth <- unique(count$Common_TRB_depth)
  density <- data.frame(Group = factor(count$Group, levels = c("Vehicle", "SLC-3010")),
                        Label = count$Sample_label, x = count$CD8_Tcells_million_per_g,
                        y = count$AH1_expected_richness_at_common_depth)
  # Panel e compares observed richness per mouse; panel f uses the rarefied expectation.
  group_levels <- c("Vehicle", "SLC-3010")
  group_colors <- c("Vehicle" = "#6E6E6E", "SLC-3010" = "#D94141")
  count$Group <- factor(count$Group, levels = group_levels)
  count_plot <- ggplot2::ggplot(count, ggplot2::aes(Group, AH1_Exact_CDR3_Matches_Any, colour = Group)) +
      ggplot2::geom_boxplot(ggplot2::aes(fill = Group), alpha = 0.13, width = 0.48, outlier.shape = NA) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.055, height = 0, seed = 3010), size = 2.5) +
      ggplot2::scale_colour_manual(values = group_colors, guide = "none") +
      ggplot2::scale_fill_manual(values = group_colors, guide = "none") +
      ggplot2::expand_limits(y = 0) + ggplot2::theme_classic(base_size = 10) +
      ggplot2::labs(x = NULL, y = "Detected AH1-reference CDR3\u03b2 sequences (n)",
          caption = sprintf("Each point is a tumor. Unpaired two-sided Student's t-test P = %s.",
                            format.pval(result$tables$student_comparison$p, digits = 3))) +
      ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0))
  density_plot <- ggplot2::ggplot(density, ggplot2::aes(x, y, colour = Group)) +
      ggplot2::geom_point(size = 2.7) +
      ggrepel::geom_text_repel(ggplot2::aes(label = Label), seed = 3010, size = 3, max.overlaps = Inf,
                               min.segment.length = 0, show.legend = FALSE) +
      ggplot2::scale_colour_manual("Treatment", values = group_colors, breaks = group_levels) +
      ggplot2::expand_limits(x = 0, y = 0) + ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(legend.position = "bottom", plot.caption = ggplot2::element_text(hjust = 0)) +
      ggplot2::labs(x = "Intratumoral CD8 T cells (\u00d710\u2076/g tumor)",
          y = "Expected AH1-reference CDR3\u03b2 sequences (n)",
          caption = sprintf("Richness is evaluated at a common depth of %s eligible TRB counts.", format(common_depth)))
  frequencies <- result$tables$sequence_frequencies
  heat <- frequencies[frequencies$CDR3_aa %in% result$detected_cdr3, , drop = FALSE]
  if (!nrow(heat)) stop("No detected reference sequences to display in the heatmap.")
  positive <- heat$Heat_TRB[heat$Heat_TRB > 0]
  cap <- as.numeric(stats::quantile(positive, 0.95, type = 7, names = FALSE))
  heat$Sample_label <- factor(heat$Sample_label, levels = count$Sample_label)
  heat$CDR3_aa <- factor(heat$CDR3_aa, levels = rev(result$detected_cdr3))
  heat$Group <- factor(heat$Group, levels = group_levels)
  heatmap <- ggplot2::ggplot(heat, ggplot2::aes(Sample_label, CDR3_aa, fill = Heat_TRB)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.15) +
      ggplot2::facet_grid(cols = ggplot2::vars(Group), scales = "free_x", space = "free_x") +
      ggplot2::scale_fill_gradientn("log10(TRB frequency per million + 1)",
          colours = c("#F2F2F2", "#440154", "#3B528B", "#21918C", "#5EC962", "#FDE725"),
          values = c(0, .002, .20, .45, .72, 1), limits = c(0, cap), oob = scales::squish,
          breaks = c(0, cap / 2, cap), labels = c("0", sprintf("%.2f", cap / 2), sprintf("\u2265%.2f", cap))) +
      ggplot2::labs(x = "Tumor", y = "Matched AH1-reference CDR3\u03b2",
          caption = "Grey denotes zero. Colour saturates at the 95th percentile of positive values; the sequence table gives the full frequencies.") +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "bottom",
                     axis.text.y = ggplot2::element_text(size = 7), plot.caption = ggplot2::element_text(hjust = 0))
  list("Figure 3e | Observed AH1 richness" = count_plot,
       "Figure 3f | CD8 density and depth-adjusted richness" = density_plot,
       "Figure 3g | AH1-matched clonotypes" = heatmap)
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
  cat("Fig. 3e-g | AH1-reference TCR sequences detected per tumor\nUnpaired, two-sided Student t-test with equal variances on observed richness; no multiplicity adjustment. Common-depth richness is a separate descriptive estimate.\n")
  show("student_comparison", "Observed richness: Student t-test and 95% confidence interval", c("group1", "group2", "n1", "n2", "mean1", "mean2", "vehicle_minus_slc", "ci95_low", "ci95_high", "t", "df", "p"))
  show("per_tumor", "Observed and common-depth richness per tumor", c("Sample_label", "Group", "TRB_Count", "AH1_richness", "Common_TRB_depth", "AH1_expected_richness_at_common_depth"))
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
    paths <- c(frames[basename(frames) == "Fig_03E-G_AH1_reference_TCR.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_03E-G_AH1_reference_TCR.R"], join(getwd(), "Fig_03E-G_AH1_reference_TCR.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_03E-G_AH1_reference_TCR.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply one workbook path or directory, or use Source in RStudio.", call. = FALSE)
    path <- path.expand(if (length(args)) args[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR"))
    book <- "Source_Data_TCR_inputs.xlsx"
    sheets <- c("AH1_reference", "S1_clone_counts", "S2_clone_counts", "S3_clone_counts", "Sample_information", "V1_clone_counts", "V2_clone_counts", "V3_clone_counts", "V4_clone_counts")
    expected <- "1d2bfb68a62e593cc2500217fc1d73b2df839afbabd3b35911b093568ab9509f"
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
      run <- join(base, "Fig_03E-G_AH1_reference_TCR", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
      needed <- c("ggplot2", "ggrepel", "scales")
      missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
      if (length(missing)) stop("Install these plotting packages: ", paste(missing, collapse = ", "))
      if (interactive() && identical(Sys.getenv("RSTUDIO"), "1") && !identical(names(grDevices::dev.cur()), "RStudioGD"))
        grDevices::dev.new(noRStudioGD = FALSE)
      measure <- !interactive() && grDevices::dev.cur() == 1L
      if (measure) grDevices::pdf(NULL)
      plots <- tryCatch(suppressPackageStartupMessages(plot_results(analysis_result)), finally = if (measure) grDevices::dev.off())
      for (i in seq_along(plots)) {
        size <- analysis_result$sizes[[names(plots)[i]]]
        if (is.null(size)) size <- c(9, 7)
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
