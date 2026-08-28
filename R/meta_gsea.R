POSITIVE_COLOUR <- "#077F97"
NEGATIVE_COLOUR <- "#400257"

read_dge_cohort <- function(path) {
  cohort <- tools::file_path_sans_ext(basename(path))
  data <- data.table::fread(path, data.table = FALSE)
  required <- c("gene_id", "log2FoldChange")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(cohort, " is missing required column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(data$gene_id)) {
    stop(cohort, " contains duplicated gene_id values. Run handleGeneDups ",
         "before iGSEA.", call. = FALSE)
  }

  if (!"lfcSE" %in% names(data)) {
    if (!"stat" %in% names(data)) {
      stop(cohort, " requires lfcSE, or stat to derive lfcSE.", call. = FALSE)
    }
    data$lfcSE <- ifelse(
      is.finite(data$stat) & data$stat != 0,
      abs(data$log2FoldChange / data$stat), NA_real_)
  }

  data$gene_id <- as.character(data$gene_id)
  data$log2FoldChange <- suppressWarnings(as.numeric(data$log2FoldChange))
  data$lfcSE <- suppressWarnings(as.numeric(data$lfcSE))
  data <- data[nzchar(data$gene_id) & is.finite(data$log2FoldChange) &
                 is.finite(data$lfcSE) & data$lfcSE > 0,
               c("gene_id", "log2FoldChange", "lfcSE"), drop = FALSE]
  data$cohort <- cohort
  data
}

read_dge_directory <- function(input_dir) {
  files <- sort(list.files(input_dir, pattern = "\\.csv$", full.names = TRUE,
                           ignore.case = TRUE))
  if (length(files) < 2L) {
    stop("iGSEA requires at least two cohort CSV files in ", input_dir,
         call. = FALSE)
  }
  cohorts <- lapply(files, read_dge_cohort)
  names(cohorts) <- tools::file_path_sans_ext(basename(files))
  dplyr::bind_rows(cohorts)
}

fixed_effect <- function(beta, se) {
  weights <- 1 / se^2
  estimate <- sum(weights * beta) / sum(weights)
  estimate_se <- sqrt(1 / sum(weights))
  z <- estimate / estimate_se
  c(beta_fixed = estimate, se_fixed = estimate_se, z_fixed = z,
    p_fixed = 2 * stats::pnorm(-abs(z)))
}

random_effect <- function(beta, se, method) {
  if (length(beta) < 2L) {
    return(c(beta_random = NA_real_, se_random = NA_real_,
             z_random = NA_real_, p_random = NA_real_, tau2 = NA_real_,
             Q = NA_real_, Q_pvalue = NA_real_, I2 = NA_real_))
  }
  fit <- tryCatch(
    metafor::rma.uni(yi = beta, sei = se, method = method, test = "z"),
    error = function(error) NULL
  )
  if (is.null(fit)) {
    return(c(beta_random = NA_real_, se_random = NA_real_,
             z_random = NA_real_, p_random = NA_real_, tau2 = NA_real_,
             Q = NA_real_, Q_pvalue = NA_real_, I2 = NA_real_))
  }
  c(beta_random = unname(fit$b[1]), se_random = fit$se,
    z_random = fit$zval, p_random = fit$pval, tau2 = fit$tau2,
    Q = fit$QE, Q_pvalue = fit$QEp, I2 = fit$I2)
}

meta_analyse_gene <- function(data, method) {
  fixed <- fixed_effect(data$log2FoldChange, data$lfcSE)
  random <- random_effect(data$log2FoldChange, data$lfcSE, method)
  effects <- data$log2FoldChange[data$log2FoldChange != 0]
  agreement <- if (length(effects)) {
    max(mean(effects > 0), mean(effects < 0))
  } else {
    NA_real_
  }
  data.frame(
    gene_id = data$gene_id[[1]], k = nrow(data),
    as.list(fixed), as.list(random), direction_agreement = agreement,
    check.names = FALSE
  )
}

run_gene_meta_analysis <- function(long_data, method = "REML") {
  method <- toupper(method)
  supported <- c("REML", "ML", "DL", "HE", "HS", "SJ", "EB", "PM")
  if (!method %in% supported) {
    stop("Unsupported heterogeneity estimator: ", method,
         ". Choose one of: ", paste(supported, collapse = ", "), call. = FALSE)
  }
  groups <- split(long_data, long_data$gene_id)
  result <- dplyr::bind_rows(lapply(groups, meta_analyse_gene, method = method))
  result$padj_fixed <- stats::p.adjust(result$p_fixed, method = "BH")
  result$padj_random <- stats::p.adjust(result$p_random, method = "BH")
  result[order(result$gene_id), , drop = FALSE]
}

make_ranked_vector <- function(meta_results, model = "random", min_cohorts = 2L) {
  statistic <- switch(
    model,
    random = "z_random",
    fixed = "z_fixed",
    stop("model must be 'fixed' or 'random'", call. = FALSE)
  )
  keep <- meta_results$k >= min_cohorts & is.finite(meta_results[[statistic]])
  ranked <- meta_results[keep, c("gene_id", statistic), drop = FALSE]
  values <- ranked[[statistic]]
  names(values) <- ranked$gene_id
  values <- sort(values, decreasing = TRUE)
  if (anyDuplicated(values)) {
    values <- values + (seq_along(values) %% 1000L) * 1e-12
  }
  values
}

run_fgsea_meta <- function(ranks, gmt_path, padj_threshold = 0.05,
                           min_size = 15L, max_size = 2000L) {
  pathways <- fgsea::gmtPathways(gmt_path)
  result <- fgsea::fgseaMultilevel(
    pathways = pathways, stats = ranks,
    minSize = min_size, maxSize = max_size, eps = 0,
    nproc = 1L,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
  )
  result <- as.data.frame(result)
  result <- result[order(result$padj), , drop = FALSE]
  result$Direction <- ifelse(
    result$NES > 0, "Upregulated", "Downregulated")
  if ("leadingEdge" %in% names(result)) {
    result$leadingEdge <- vapply(
      result$leadingEdge, paste, collapse = ", ", FUN.VALUE = character(1))
  }
  list(all = result,
       significant = result[is.finite(result$padj) &
                              result$padj < padj_threshold, , drop = FALSE])
}

plot_meta_gsea <- function(result, contrast) {
  ggplot2::ggplot(
    result,
    ggplot2::aes(x = stats::reorder(pathway, NES), y = NES,
                 fill = Direction)
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Pathway", y = "Normalized Enrichment Score",
                  title = contrast) +
    ggplot2::theme(
      text = ggplot2::element_text(size = 8),
      title = ggplot2::element_text(size = 18),
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::scale_fill_manual(values = c(
      "Upregulated" = POSITIVE_COLOUR,
      "Downregulated" = NEGATIVE_COLOUR
    )) +
    ggprism::theme_prism(base_size = 12)
}

write_igsea_outputs <- function(meta_results, fgsea_results, plot,
                                output_dir, contrast) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(
    meta_results, file.path(output_dir, "meta_gene_level_results.csv"))
  data.table::fwrite(
    fgsea_results$significant,
    file.path(output_dir, "meta_fgsea_results.csv"))
  data.table::fwrite(
    fgsea_results$all,
    file.path(output_dir, "meta_fgsea_results_all.csv"))

  grDevices::png(file.path(output_dir, paste0(contrast, ".png")),
                 width = 3200, height = 1800, res = 300)
  print(plot)
  grDevices::dev.off()
  grDevices::pdf(file.path(output_dir, paste0(contrast, ".pdf")),
                 width = 12, height = 10)
  print(plot)
  grDevices::dev.off()
}

run_igsea <- function(input_dir, output_dir, gmt_path, contrast,
                      model = "random", method = "REML", min_cohorts = 2L,
                      padj_threshold = 0.05, min_size = 15L,
                      max_size = 2000L) {
  long <- read_dge_directory(input_dir)
  meta <- run_gene_meta_analysis(long, method = method)
  ranks <- make_ranked_vector(meta, model = model, min_cohorts = min_cohorts)
  if (length(ranks) < min_size) {
    stop("Only ", length(ranks), " genes passed the minimum cohort rule; ",
         "not enough for GSEA min-size ", min_size, ".", call. = FALSE)
  }
  enrichment <- run_fgsea_meta(
    ranks, gmt_path, padj_threshold, min_size, max_size)
  plot_data <- enrichment$significant
  if (!nrow(plot_data)) plot_data <- enrichment$all[0, , drop = FALSE]
  plot <- plot_meta_gsea(plot_data, contrast)
  write_igsea_outputs(meta, enrichment, plot, output_dir, contrast)
  invisible(list(meta = meta, enrichment = enrichment, ranks = ranks))
}
