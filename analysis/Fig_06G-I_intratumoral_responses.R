# Fig. 6g-i | Intratumoral immune-cell comparisons
# Compare mouse measurements with four two-sided, unpaired equal-variance t-tests per panel
# and apply Holm adjustment within each panel.
# Panels G and H report millions of cells per gram of tumor. Panel I reports the MP CD8 T-cell
# to Treg ratio.
#
# Input: put Source_Data.xlsx beside this script or in the repository's data/ folder.
# Numerical packages: digest, readxl.
# Plotting packages: ggbeeswarm, ggplot2, patchwork.
# In RStudio, open this file and click Source to save the CSVs and then display the plots.
# To see the result, enter these lines in the RStudio Console after Source:
#   analysis_result$tables$comparisons
#   analysis_result$output_dir   # Folder containing Fig_06G-I_comparisons.csv
# For tables only, enter options(SLC3010.tables_only = TRUE) in the Console, then click Source.
# To include plots again, enter options(SLC3010.tables_only = NULL) and click Source.

# Figure 6g–i. Intratumoral cell counts (1e6 cells per gram of tumour) and MP CD8 T-cell/Treg ratio.
# The two-sided unpaired Student's t-tests use the original, untransformed values and equal variances.
# Each panel reports four comparisons with Holm-adjusted P values. The symbols are calculated from
# those adjusted values. The g and i axes are logarithmic. Whiskers use 1.5 times the interquartile range.

# This stage calculates and saves the numerical tables before returning them.
analyze_data <- function(workbook, output_dir) {
  # Excel ranges used for named datasets in this analysis.
  table_ranges <- c(
    "Fig_06G_intratumoral_MP_CD8" = "'Figure 6'!A230:E268",
    "Fig_06H_intratumoral_Treg" = "'Figure 6'!A274:E312",
    "Fig_06I_intratumoral_MP_CD8_Treg_ratio" = "'Figure 6'!A318:E356"
  )

  dataset_ids <- c("Fig_06G_intratumoral_MP_CD8", "Fig_06H_intratumoral_Treg",
                   "Fig_06I_intratumoral_MP_CD8_Treg_ratio")
  raw <- do.call(rbind, lapply(dataset_ids, function(id) {
    as.data.frame(readxl::read_excel(workbook, range = table_ranges[[id]], col_types = c("numeric", "text", "text", "text", "numeric"), .name_repair = "minimal"))
  }))

  needed <- c("mouse_id", "strain", "treatment", "readout", "value")
  if (!all(needed %in% names(raw))) stop("Missing columns: ", paste(setdiff(needed,names(raw)),collapse=", "))
  raw <- raw[needed]
  if (anyNA(raw)) stop("Missing value/identifier detected; exclusions must be reviewed, not silently removed.")
  if (!is.numeric(raw$value) || any(!is.finite(raw$value))) stop("value must be finite numeric data.")
  raw$mouse_id <- as.character(raw$mouse_id)
  raw$strain <- as.character(raw$strain)
  raw$treatment <- as.character(raw$treatment)
  raw$readout <- as.character(raw$readout)
  groups <- c("wt_PBS", "wt_SLC-3010", "triple hIL-2_PBS", "triple hIL-2_SLC-3010")
  readouts <- c("MP_CD8_T", "Treg", "MP_CD8_T/Treg ratio")
  panel_letters <- setNames(c("g","h","i"), readouts)
  raw$group_id <- paste(raw$strain,raw$treatment,sep="_")
  if (any(!raw$group_id %in% groups) || any(!raw$readout %in% readouts)) stop("Unknown group/readout.")
  if (anyDuplicated(raw[c("mouse_id","strain","treatment","readout")])) stop("Duplicate mouse/readout.")

  counts <- table(factor(raw$readout,levels=readouts),factor(raw$group_id,levels=groups))
  if (any(counts < 2)) stop("At least two observations per group/readout are required.")

  raw$strain <- factor(raw$strain, levels=c("wt","triple hIL-2"))
  raw$treatment <- factor(raw$treatment,levels=c("PBS","SLC-3010"))
  raw$group_id <- factor(raw$group_id,levels=groups)

  summary_rows <- list()
  for (r in readouts) for (g in groups) {
    v <- raw$value[raw$readout==r & raw$group_id==g]
    summary_rows[[length(summary_rows)+1]] <- data.frame(panel=panel_letters[[r]],readout=r,group=g,
      n=length(v),mean=mean(v),sd=sd(v),median=median(v),min=min(v),max=max(v))
  }
  group_summary <- do.call(rbind,summary_rows)

  contrast <- data.frame(
    comparison=c("within_wt","within_triple","between_vehicle","between_slc"),
    contrast=c("A vs B","C vs D","A vs C","B vs D"),
    group1=groups[c(1,3,1,2)],group2=groups[c(2,4,3,4)],stringsAsFactors=FALSE)
  rows <- list()
  for (r in readouts) for (j in seq_len(nrow(contrast))) {
    cmp <- contrast[j,,drop=FALSE]
    x <- raw$value[raw$readout==r & raw$group_id==cmp$group1]
    y <- raw$value[raw$readout==r & raw$group_id==cmp$group2]
    # P-values always from two-sided, unpaired, equal-variance Student t-tests.
    tt <- stats::t.test(y, x, alternative="two.sided", paired=FALSE, var.equal=TRUE)
    rows[[length(rows)+1]] <- cbind(
      data.frame(panel=panel_letters[[r]], readout=r), cmp,
      data.frame(n1=length(x), n2=length(y), mean1=mean(x), mean2=mean(y),
        sd1=sd(x), sd2=sd(y), mean_difference=mean(y)-mean(x),
        fold_mean2_over_mean1=mean(y)/mean(x),
        student_t=unname(tt$statistic), student_df=unname(tt$parameter),
        student_p=tt$p.value, student_ci_low=tt$conf.int[1],
        student_ci_high=tt$conf.int[2])
    )
  }
  results <- do.call(rbind, rows)
  stopifnot(nrow(results) == 12L, all(table(results$panel) == 4L))
  # Four comparisons form a correction family separately for g, h and i.
  results$p_display <- ave(results$student_p, results$panel,
                           FUN=function(z) stats::p.adjust(z, method="holm"))
  p_to_star <- function(p) {
    ifelse(p < 0.0001, "****", ifelse(p < 0.001, "***",
           ifelse(p < 0.005, "**", ifelse(p < 0.05, "*", "ns"))))
  }
  results$star <- p_to_star(results$p_display)
  results$test <- "two-sided unpaired Student t-test"
  results$p_adjustment <- "Holm within panel (four comparisons)"

  utils::write.csv(raw, file.path(output_dir, "Fig_06G-I_individual_data.csv"), row.names = FALSE)
  group_metadata <- unique(raw[c("group_id", "strain", "treatment")])
  group_summary$group_id <- group_summary$group
  group_summary$strain <- group_metadata$strain[match(group_summary$group, group_metadata$group_id)]
  group_summary$treatment <- group_metadata$treatment[match(group_summary$group, group_metadata$group_id)]
  summary_columns <- c("readout", "strain", "treatment", "group_id", "n", "mean", "sd", "median", "min", "max")
  utils::write.csv(group_summary[summary_columns], file.path(output_dir, "Fig_06G-I_summary.csv"), row.names = FALSE)
  utils::write.csv(results, file.path(output_dir, "Fig_06G-I_comparisons.csv"), row.names = FALSE)

  list(tables = list(comparisons = results, summary = group_summary[summary_columns], observations = raw))
}

# Draw the panels from the results returned by analyze_data().
plot_results <- function(result) {
  raw <- result$tables$observations
  results <- result$tables$comparisons
  groups <- levels(raw$group_id)
  readouts <- unique(results$readout)
  panel_letters <- setNames(c("g", "h", "i"), readouts)
  if (any(raw$value[raw$readout %in% readouts[c(1, 3)]] <= 0))
    stop("The logarithmic panels require positive values for display; numerical results are available.")

  BASE_FAMILY <- "Arial"
  BASE_SIZE_PT <- 5.75
  Y_TITLE_SIZE_PT <- 5.80
  X_TEXT_SIZE_PT <- 5.50
  Y_TEXT_SIZE_PT <- 5.50
  STAT_TEXT_SIZE_MM <- 1.75
  LEGEND_TEXT_SIZE_PT <- 6.00
  POINT_SIZE_MM <- 1.25
  AXIS_LINE_WIDTH <- 0.34
  BOX_LINE_WIDTH <- 0.38
  # Bracket stroke, tip length and label spacing in millimetres.
  BRACKET_LINE_WIDTH_MM <- 0.24
  BRACKET_TIP_MM <- 0.38
  BRACKET_TEXT_GAP_MM <- 0.40
  BOXPLOT_COEF <- 1.5  # ggplot2 default: whiskers extend to the most extreme value within 1.5 x IQR

  # Bracket levels are fractions of panel height on both linear and log axes.
  # Within-strain comparisons share a level below the two cross-strain comparisons.
  BRACKET_LANES <- c(within_wt = 0.745, within_triple = 0.745,
                     between_vehicle = 0.835, between_slc = 0.925)
  # Lower bounds and minimum upper bounds for each readout.
  Y_MIN <- c(MP_CD8_T = 0.1, Treg = 0, "MP_CD8_T/Treg ratio" = 0.1)
  Y_UPPER_MIN <- c(MP_CD8_T = 600, Treg = 5, "MP_CD8_T/Treg ratio" = 1500)
  # Reserve the upper part of the panel for the comparison brackets.
  MAX_DATA_NPC <- 0.69

  REQUESTED_FONT <- BASE_FAMILY
  if (requireNamespace("systemfonts", quietly=TRUE)) {
    font_names <- unique(systemfonts::system_fonts()$family)
    if (!BASE_FAMILY %in% font_names) {
      if ("Nimbus Sans" %in% font_names) BASE_FAMILY <- "Nimbus Sans"
      else if ("Liberation Sans" %in% font_names) BASE_FAMILY <- "Liberation Sans"
      else BASE_FAMILY <- "sans"
      warning("Font ", REQUESTED_FONT, " unavailable; using ", BASE_FAMILY, ".", call.=FALSE)
    }
  }

  x_positions <- setNames(c(1,2,3.3,4.3),groups)
  raw$x_num <- unname(x_positions[raw$group_id])

  colors <- c(PBS="#111111", "SLC-3010"="#FF1A1A")
  shapes <- c(wt=16, "triple hIL-2"=17)
  y_labels <- c(MP_CD8_T="MP CD8 T (×10⁶ cells/g tumor)",Treg="Treg (×10⁶ cells/g tumor)",
                "MP_CD8_T/Treg ratio"="MP CD8 T / Treg ratio")
  PT_PER_MM <- 72.27/25.4
  X_LIMITS <- c(0.45,4.85)

  # Brackets are a grid overlay in DISPLAY coordinates. Each is a four-point polyline:
  # tip -> left corner -> right corner -> tip. Tips and label gaps use physical mm.
  make_bracket_grob <- function(s) {
    grobs <- list()
    for (k in seq_len(nrow(s))) {
      a <- s[k,,drop=FALSE]
      x1 <- (x_positions[[a$group1]]-X_LIMITS[1])/diff(X_LIMITS)
      x2 <- (x_positions[[a$group2]]-X_LIMITS[1])/diff(X_LIMITS)
      y <- grid::unit(BRACKET_LANES[[a$comparison]],"npc")
      tip <- grid::unit(BRACKET_TIP_MM,"mm")
      grobs[[length(grobs)+1]] <- grid::polylineGrob(
        x=grid::unit(c(x1,x1,x2,x2),"npc"), y=grid::unit.c(y-tip,y,y,y-tip),
        gp=grid::gpar(col="black",fill=NA,lwd=BRACKET_LINE_WIDTH_MM*PT_PER_MM,
                      lineend="butt",linejoin="mitre"))
      grobs[[length(grobs)+1]] <- grid::textGrob(a$star,
        x=grid::unit((x1+x2)/2,"npc"),y=y+grid::unit(BRACKET_TEXT_GAP_MM,"mm"),
        just=c("centre","bottom"),gp=grid::gpar(fontfamily=BASE_FAMILY,
          fontsize=STAT_TEXT_SIZE_MM*PT_PER_MM,fontface="bold",col="black"))
    }
    do.call(grid::grobTree,grobs)
  }

  base_theme <- function() {
    ggplot2::theme_classic(base_size=BASE_SIZE_PT,base_family=BASE_FAMILY)+
      ggplot2::theme(plot.title=ggplot2::element_blank(),plot.tag=ggplot2::element_blank(),
        axis.title.x=ggplot2::element_blank(),
        axis.title.y=ggplot2::element_text(size=Y_TITLE_SIZE_PT,face="plain",color="black",
                                         margin=ggplot2::margin(r=1.2)),
        axis.text.x=ggplot2::element_text(size=X_TEXT_SIZE_PT,color="black",angle=90,
                   hjust=1,vjust=0.5,margin=ggplot2::margin(t=1.0)),
        axis.text.y=ggplot2::element_text(size=Y_TEXT_SIZE_PT,color="black"),
        axis.line=ggplot2::element_line(linewidth=AXIS_LINE_WIDTH,color="black"),
        axis.ticks=ggplot2::element_line(linewidth=AXIS_LINE_WIDTH,color="black"),
        axis.ticks.length=grid::unit(1.6,"pt"),legend.position="none",
        plot.margin=ggplot2::margin(t=3,r=1,b=1,l=1,unit="pt"))
  }
  make_panel <- function(r,s) {
    dat <- raw[raw$readout==r,,drop=FALSE]
    log_panel <- r %in% c("MP_CD8_T","MP_CD8_T/Treg ratio")
    lower <- Y_MIN[[r]]
    # Scale the upper limit to leave room between the largest observation and the lowest bracket.
    needed_upper <- if (log_panel) 10^(log10(lower)+(log10(max(dat$value))-log10(lower))/MAX_DATA_NPC)
                    else lower+(max(dat$value)-lower)/MAX_DATA_NPC
    upper <- max(Y_UPPER_MIN[[r]],needed_upper)
    if (min(dat$value)<lower) stop("Data below configured lower limit; revise Y_MIN explicitly.")
    p <- ggplot2::ggplot(dat,ggplot2::aes(x=x_num,y=value))+
      ggplot2::stat_boxplot(ggplot2::aes(group=group_id,color=treatment),geom="errorbar",
                            width=0.34,coef=BOXPLOT_COEF,linewidth=AXIS_LINE_WIDTH)+
      ggplot2::geom_boxplot(ggplot2::aes(group=group_id,color=treatment),width=0.52,
                            fill="white",coef=BOXPLOT_COEF,outlier.shape=NA,linewidth=BOX_LINE_WIDTH)
    p <- p+ggbeeswarm::geom_quasirandom(
        ggplot2::aes(color=treatment,shape=strain,group=group_id),width=0.12,
        size=POINT_SIZE_MM,stroke=0.35,alpha=1)
    p <- p+ggplot2::scale_x_continuous(breaks=unname(x_positions),
            labels=c("PBS","SLC-3010","PBS","SLC-3010"),limits=X_LIMITS,
            expand=ggplot2::expansion(mult=0))+
      ggplot2::scale_color_manual(values=colors,drop=FALSE)+
      ggplot2::scale_shape_manual(values=shapes,drop=FALSE)+
      ggplot2::labs(y=y_labels[[r]])+base_theme()
    if (log_panel) {
      br <- c(0.1,1,10,100,1000,10000); br <- br[br>=lower & br<=upper]
      p <- p+ggplot2::scale_y_log10(breaks=br,labels=vapply(br,function(z) format(z,scientific=FALSE,trim=TRUE),character(1)),
                                   expand=ggplot2::expansion(mult=0))+
        ggplot2::annotation_logticks(sides="l",short=grid::unit(0.55,"pt"),
          mid=grid::unit(0.90,"pt"),long=grid::unit(1.60,"pt"),color="black",linewidth=0.22)
    } else p <- p+ggplot2::scale_y_continuous(breaks=seq(0,floor(upper)),
                                            expand=ggplot2::expansion(mult=0))
    p+ggplot2::coord_cartesian(ylim=c(lower,upper),expand=FALSE,clip="off")+
      ggplot2::annotation_custom(make_bracket_grob(s[s$readout==r,,drop=FALSE]),
        xmin=X_LIMITS[1],xmax=X_LIMITS[2],ymin=lower,ymax=upper)
  }
  make_strain_legend <- function() {
    ggplot2::ggplot(data.frame(x=c(0.23,0.59),y=0.5,strain=c("wt","triple hIL-2"),
                               label=c("C57BL/6","Triple hIL-2")))+
      ggplot2::geom_point(ggplot2::aes(x=x,y=y,shape=strain),size=1.6)+
      ggplot2::geom_text(ggplot2::aes(x=x+0.035,y=y,label=label),hjust=0,
                         size=LEGEND_TEXT_SIZE_PT/PT_PER_MM,family=BASE_FAMILY)+
      ggplot2::scale_shape_manual(values=shapes)+
      ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand=FALSE,clip="off")+
      ggplot2::theme_void()+ggplot2::theme(legend.position="none",
        plot.margin=ggplot2::margin(t=2.5,r=0,b=0,l=0,unit="pt"))
  }
  make_combined <- function(s) {
    panels <- lapply(readouts,make_panel,s=s)
    row <- patchwork::wrap_plots(panels,nrow=1,widths=c(1,1,1))
    patchwork::wrap_plots(row,make_strain_legend(),ncol=1,heights=c(0.90,0.10))
  }

  panels <- make_combined(results)

  list("Fig. 6g–i | Intratumoral responses" = panels)
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
  cat("Fig. 6g-i | Intratumoral immune-cell comparisons\nFour unpaired, two-sided Student t-tests per panel; Holm adjustment is within each panel. p_display is adjusted; confidence intervals are unadjusted. Mean differences and t statistics are group 2 minus group 1.\n")
  show("comparisons", "Mean differences and unadjusted 95% confidence intervals", c("panel", "readout", "group1", "group2", "n1", "n2", "mean_difference", "student_ci_low", "student_ci_high"))
  show("comparisons", "Student t-tests and Holm-adjusted P values", c("panel", "group1", "group2", "student_t", "student_df", "student_p", "p_display"))
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
    paths <- c(frames[basename(frames) == "Fig_06G-I_intratumoral_responses.R"], if (direct) cli, editor[basename(as.character(editor)) == "Fig_06G-I_intratumoral_responses.R"], join(getwd(), "Fig_06G-I_intratumoral_responses.R"), frames[grepl("[.]R$", frames, ignore.case = TRUE)])
    paths <- paths[!is.na(paths) & file.exists(paths) & !dir.exists(paths)]
    here <- if (length(paths)) dirname(normalizePath(paths[1L], winslash = "/")) else getwd()
    parents <- c(dirname(here), here)
    repo <- head(parents[file.exists(join(parents, "README.md")) & file.exists(join(parents, "analysis", "Fig_06G-I_intratumoral_responses.R"))], 1L)
    args <- if (direct) commandArgs(trailingOnly = TRUE) else character()
    tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
    if (length(args) > 1L || any(startsWith(args, "--"))) stop("Supply one workbook path or directory, or use Source in RStudio.", call. = FALSE)
    path <- path.expand(if (length(args)) args[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR"))
    book <- "Source_Data.xlsx"
    sheets <- c("Figure 6")
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
      run <- join(base, "Fig_06G-I_intratumoral_responses", basename(tempfile(format(Sys.time(), "%Y%m%d-%H%M%S-"))))
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
      needed <- c("ggbeeswarm", "ggplot2", "patchwork")
      missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
      if (length(missing)) stop("Install these plotting packages: ", paste(missing, collapse = ", "))
      if (interactive() && identical(Sys.getenv("RSTUDIO"), "1") && !identical(names(grDevices::dev.cur()), "RStudioGD"))
        grDevices::dev.new(noRStudioGD = FALSE)
      measure <- !interactive() && grDevices::dev.cur() == 1L
      if (measure) grDevices::pdf(NULL)
      plots <- tryCatch(suppressPackageStartupMessages(plot_results(analysis_result)), finally = if (measure) grDevices::dev.off())
      for (i in seq_along(plots)) {
        size <- analysis_result$sizes[[names(plots)[i]]]
        if (is.null(size)) size <- c(12, 5)
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
