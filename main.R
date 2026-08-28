#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  getwd()
}
source(file.path(script_dir, "R", "meta_gsea.R"))

usage <- function(status = 0L) {
  cat(
    "iGSEA - gene-level meta-analysis followed by preranked GSEA\n\n",
    "Usage:\n",
    "  iGSEA -d DGE_DIRECTORY -g GENE_SETS.gmt -c CONTRAST [options]\n\n",
    "Required:\n",
    "  -d, --data DIR          Directory containing cohort DGE CSV files\n",
    "  -g, --gene-set FILE     GMT gene-set file\n",
    "  -c, --contrast TEXT     Contrast label used as the plot title\n\n",
    "Options:\n",
    "  -o, --output DIR        Output directory [results]\n",
    "  -m, --model NAME        Ranking model: random or fixed [random]\n",
    "      --method NAME       Between-cohort variance estimator [REML]\n",
    "      --min-cohorts INT   Minimum contributing cohorts per gene [2]\n",
    "      --padj NUMBER       Significant pathway threshold [0.05]\n",
    "      --min-size INT      Minimum gene-set size [15]\n",
    "      --max-size INT      Maximum gene-set size [2000]\n",
    "  -h, --help              Show this help\n",
    "  -V, --version           Show version\n\n",
    "Methods: REML, ML, DL, HE, HS, SJ, EB, PM\n",
    sep = ""
  )
  quit(status = status)
}

parse_arguments <- function(args) {
  options <- list(
    data = NULL, gene_set = NULL, contrast = NULL, output = "results",
    model = "random", method = "REML", min_cohorts = 2L,
    padj = 0.05, min_size = 15L, max_size = 2000L
  )
  aliases <- c(
    "-d" = "data", "--data" = "data",
    "-g" = "gene_set", "--gene-set" = "gene_set",
    "-c" = "contrast", "--contrast" = "contrast",
    "-o" = "output", "--output" = "output",
    "-m" = "model", "--model" = "model",
    "--method" = "method", "--min-cohorts" = "min_cohorts",
    "--padj" = "padj", "--min-size" = "min_size",
    "--max-size" = "max_size"
  )
  index <- 1L
  while (index <= length(args)) {
    argument <- args[[index]]
    if (argument %in% c("-h", "--help")) usage()
    if (argument %in% c("-V", "--version")) {
      cat("iGSEA 1.0.0\n")
      quit(status = 0L)
    }
    if (!argument %in% names(aliases)) {
      stop("Unknown option: ", argument, call. = FALSE)
    }
    if (index == length(args)) stop("Missing value for ", argument, call. = FALSE)
    options[[aliases[[argument]]]] <- args[[index + 1L]]
    index <- index + 2L
  }

  missing <- names(options)[c("data", "gene_set", "contrast")][
    vapply(options[c("data", "gene_set", "contrast")], is.null, logical(1))]
  if (length(missing)) {
    stop("Missing required option(s): ", paste(missing, collapse = ", "),
         ". Run iGSEA --help.", call. = FALSE)
  }
  if (!dir.exists(options$data)) stop("Data directory does not exist: ", options$data)
  if (!file.exists(options$gene_set)) stop("GMT file does not exist: ", options$gene_set)
  options$model <- tolower(options$model)
  if (!options$model %in% c("fixed", "random")) {
    stop("--model must be fixed or random", call. = FALSE)
  }
  options$method <- toupper(options$method)
  options$min_cohorts <- as.integer(options$min_cohorts)
  options$padj <- as.numeric(options$padj)
  options$min_size <- as.integer(options$min_size)
  options$max_size <- as.integer(options$max_size)
  if (is.na(options$min_cohorts) || options$min_cohorts < 1L) {
    stop("--min-cohorts must be a positive integer", call. = FALSE)
  }
  if (!is.finite(options$padj) || options$padj <= 0 || options$padj > 1) {
    stop("--padj must be in (0, 1]", call. = FALSE)
  }
  if (options$min_size < 1L || options$max_size < options$min_size) {
    stop("Invalid gene-set size limits", call. = FALSE)
  }
  options
}

tryCatch({
  options <- parse_arguments(commandArgs(trailingOnly = TRUE))
  message("Reading independent cohort DGE results...")
  message("Meta-analysis: model=", options$model,
          ", heterogeneity method=", options$method)
  run_igsea(
    input_dir = options$data,
    output_dir = options$output,
    gmt_path = options$gene_set,
    contrast = options$contrast,
    model = options$model,
    method = options$method,
    min_cohorts = options$min_cohorts,
    padj_threshold = options$padj,
    min_size = options$min_size,
    max_size = options$max_size
  )
  message("iGSEA completed. Results: ", normalizePath(
    options$output, mustWork = FALSE))
}, error = function(error) {
  message("Error: ", conditionMessage(error))
  quit(status = 1L)
})
