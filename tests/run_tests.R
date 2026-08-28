#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
test_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
root <- dirname(test_dir)
source(file.path(root, "R", "meta_gsea.R"))

assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

long <- data.frame(
  gene_id = rep(c("geneA", "geneB", "geneC"), each = 3),
  log2FoldChange = c(1.0, 1.2, 0.9, -1.1, -0.8, -1.0, 1.5, -1.4, 0.1),
  lfcSE = rep(0.2, 9),
  cohort = rep(c("A", "B", "C"), 3)
)
meta <- run_gene_meta_analysis(long, method = "REML")
assert(nrow(meta) == 3L, "Expected one result per gene")
assert(meta$z_random[meta$gene_id == "geneA"] > 0,
       "geneA should have a positive random-effects statistic")
assert(meta$z_random[meta$gene_id == "geneB"] < 0,
       "geneB should have a negative random-effects statistic")
assert(meta$I2[meta$gene_id == "geneC"] > meta$I2[meta$gene_id == "geneA"],
       "Discordant geneC should be more heterogeneous than geneA")

random_ranks <- make_ranked_vector(meta, model = "random", min_cohorts = 2)
fixed_ranks <- make_ranked_vector(meta, model = "fixed", min_cohorts = 2)
assert(identical(names(random_ranks)[1], "geneA"),
       "Random ranking should retain signed meta evidence")
assert(length(fixed_ranks) == 3L, "Fixed ranking lost valid genes")

temporary <- tempfile("igsea-test-")
dir.create(temporary)
duplicate_file <- file.path(temporary, "duplicate.csv")
data.table::fwrite(data.frame(
  gene_id = c("geneA", "geneA"), log2FoldChange = c(1, 2),
  lfcSE = c(0.2, 0.2)), duplicate_file)
duplicate_error <- tryCatch({
  read_dge_cohort(duplicate_file)
  FALSE
}, error = function(error) grepl("duplicated gene_id", conditionMessage(error)))
assert(duplicate_error, "Duplicate genes should fail with an explicit error")

plot_output <- file.path(temporary, "plots")
test_plot <- ggplot2::ggplot() + ggplot2::geom_blank()
write_igsea_outputs(
  meta_results = data.frame(gene_id = "geneA"),
  fgsea_results = list(
    significant = data.frame(pathway = "setA"),
    all = data.frame(pathway = "setA")
  ),
  plot = test_plot,
  output_dir = plot_output,
  contrast = "branched"
)
assert(file.exists(file.path(plot_output, "branched.png")),
       "PNG filename should be taken from the contrast")
assert(file.exists(file.path(plot_output, "branched.pdf")),
       "PDF filename should be taken from the contrast")

cat("All iGSEA tests passed.\n")
