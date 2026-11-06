# Run all 33 study analyses, each in a fresh R session. For a single panel,
# open its file in analysis/ instead; each of those scripts uses one workbook.
# For a complete run, put Source_Data.xlsx, Source_Data_RNA_inputs.xlsx and
# Source_Data_TCR_inputs.xlsx in data/ beside SLC3010.Rproj.
# In RStudio, open this file and click Source. Each child process writes its own
# CSVs and plot files and prints its output directory in the Console.
# For CSVs only, enter options(SLC3010.tables_only = TRUE) in the Console first.
# Restore plots on later runs with options(SLC3010.tables_only = NULL).
frames <- vapply(sys.frames(), function(f) {
  path <- f[["ofile"]]
  if (!is.character(path) || length(path) != 1L || is.na(path)) "" else path
}, "")
source_file <- tail(frames[basename(frames) == "RUN_ALL.R" & file.exists(frames)], 1L)
cli_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
file <- if (length(source_file)) source_file else cli_file
root <- if (length(file) && file.exists(file[1L])) dirname(normalizePath(file[1L])) else getwd()
args <- if (length(source_file)) character() else commandArgs(trailingOnly = TRUE)
tables_only <- isTRUE(getOption("SLC3010.tables_only", FALSE))
paths <- args
if (length(paths) > 1L || any(startsWith(paths, "--"))) stop("Supply one workbook directory, or use Source in RStudio.", call. = FALSE)
folder <- if (length(paths)) paths[1L] else Sys.getenv("SLC_SOURCE_DATA_DIR", unset = file.path(root, "data"))

deps <- utils::read.csv(file.path(root, "environment", "analysis_packages.csv"), stringsAsFactors = FALSE)
stage_columns <- c("Analysis_CRAN", "Analysis_Bioconductor")
if (!tables_only) stage_columns <- c(stage_columns, "Plot_CRAN", "Plot_Bioconductor")
entries <- unlist(deps[stage_columns], use.names = FALSE)
needed <- unique(trimws(unlist(strsplit(entries[!is.na(entries) & nzchar(entries)], ",", fixed = TRUE))))
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install these packages for the requested run: ", paste(missing, collapse = ", "), call. = FALSE)

input_map <- utils::read.csv(file.path(root, "docs", "analysis_inputs.csv"), stringsAsFactors = FALSE)
for (book in unique(deps$Workbook)) {
  path <- file.path(folder, book)
  if (!file.exists(path) || dir.exists(path)) stop("Missing workbook: ", path, call. = FALSE)
  sheets <- unique(input_map$Sheet[input_map$Workbook == book])
  absent <- setdiff(sheets, readxl::excel_sheets(path))
  if (length(absent)) stop("Missing worksheets in ", book, ": ", paste(absent, collapse = ", "), call. = FALSE)
}
folder <- normalizePath(folder, winslash = "/", mustWork = TRUE)
cmd <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
scripts <- file.path(root, deps$Script)
if (!all(file.exists(scripts))) stop("Use RUN_ALL.R in the extracted repository.", call. = FALSE)
run_analysis <- function(script) {
  launcher <- tempfile("SLC3010-analysis-", fileext = ".R")
  on.exit(unlink(launcher), add = TRUE)
  writeLines(c(
    sprintf("options(SLC3010.autorun = TRUE, SLC3010.tables_only = %s)", if (tables_only) "TRUE" else "FALSE"),
    sprintf("Sys.setenv(SLC_SOURCE_DATA_DIR = %s)", encodeString(folder, quote = '"')),
    sprintf("source(%s, encoding = 'UTF-8')", encodeString(normalizePath(script, winslash = "/"), quote = '"'))
  ), launcher, useBytes = TRUE)
  status <- system2(cmd, shQuote(launcher))
  if (status != 0L) stop(basename(script), " exited with status ", status, call. = FALSE)
  invisible(status)
}
for (script in scripts) {
  message("Running ", basename(script))
  run_analysis(script)
}
