# Meta-GSEA pipeline: per-gene meta-analysis then GSEA on combined stats
# Requirements: DESeq2 per-cohort results tables (data frames or files)
# Packages: data.table, dplyr, metafor, fgsea, BiocParallel (optional), msigdbr (optional)
# Install if needed:
# install.packages(c("data.table","dplyr","metafor","fgsea","BiocParallel"))
# Bioconductor packages: BiocManager::install("DESeq2"); BiocManager::install("fgsea")
pkgs <- c("data.table", "dplyr", "metafor", "fgsea", "BiocParallel", "here", "pipecore", "msigdbr", "ggplot2", "ggprism")
suppressMessages(pipecore::handleRequirements(pkgs))

# library(data.table)
# library(dplyr)
# library(metafor)
# library(fgsea)
# library(BiocParallel)  # optional for parallelization

# ------------------------
# 1) Load per-cohort DE results
# ------------------------
# Expected structure for each file / data.frame: gene, log2FoldChange, lfcSE (or stat), pvalue
# Example: a list of CSVs, one per cohort
args      <- commandArgs(TRUE)

inpath    <- normalizePath(args[1])
outpath   <- normalizePath(args[2], mustWork = F)
gmt.path  <- args[3]
contrast  <- args[4]

# print(outpath)
# inpath <- "/home/jan1/Documents/Cancer_Studies_PhD/Deciphering_Oncogenes/Workflow/PhD_Projects/Yue_BAP1/iGSEA/iGSEA_input"
# # inpath <- "/home/jan1/Documents/Cancer_Studies_PhD/Studies/Study_Phylogenetics/results/iGSEA/archive/LvB/data/"
# outpath <- "/home/jan1/Documents/Cancer_Studies_PhD/Deciphering_Oncogenes/Workflow/PhD_Projects/Yue_BAP1/iGSEA/results/"
# contrast <- "Mutated"

dir.create(outpath, F, T)
oupath <- getwd()
print(outpath)

cohort_files <- list.files(inpath, pattern="*.csv", full.names=TRUE) # adjust folder
cohort_names <- gsub("\\.csv$","", basename(cohort_files))

read_cohort <- function(file, cohort_name) {
  cat(crayon::green("Reading in ", cohort_name, "...\n"))
  df <- fread(file) %>% as.data.frame()

  # make sure column names exist; adapt if your files have different colnames
  # require 'gene' column
  if(!"gene_id" %in% colnames(df)) stop(paste("No 'gene' column in", file))
  # try to ensure needed columns
  if(!"lfcSE" %in% colnames(df) & "stat" %in% colnames(df) & "log2FoldChange" %in% colnames(df)) {
    df$lfcSE <- with(df, ifelse(stat==0 | is.na(stat), NA, log2FoldChange / stat))
  }
  if(!"lfcSE" %in% colnames(df)) {
    warning("lfcSE not found; will try to compute from log2FoldChange/stat; some SEs may be NA")
    df$lfcSE <- NA
  }
  df$cohort <- cohort_name

  df <- df %>% dplyr::select(gene_id, log2FoldChange, lfcSE, stat, pvalue, cohort)
  return(df)
}

cohort_list <- mapply(read_cohort, cohort_files, cohort_names, SIMPLIFY = FALSE)

# Check for duplicates
check_dups <- function(df){
  dups <- length(which(duplicated(df$gene_id)))
  if(dups > 0){
    df <- df[-which(duplicated(df$gene_id)),]
    return(df)
  }
  else{
    return(df)
  }
}

cohort_list <- lapply(cohort_list, check_dups)


# Combine into long table
long <- bind_rows(cohort_list)

# Optional QC: number of genes per cohort
table(long$cohort)

# ------------------------
# 2) Harmonize genes: build gene universe
# ------------------------
# Option A: If all cohorts use same gene IDs (e.g., gene symbols), proceed.
# Option B: If different (ENSEMBL vs symbol) or species differ, map orthologs first (use biomaRt or org.* packages).
# For simplicity assume same gene IDs:
cat(crayon::green("Solving gene universe...\n"))
genes_all <- sort(unique(long$gene))
length(genes_all)

dups <- long |>
  dplyr::count(gene_id, cohort, name = "n") |>
  dplyr::filter(n > 1)

cat(crayon::green("No duplicated genes were found...\n"))
dups |> dplyr::arrange(dplyr::desc(n))

# ------------------------
# 3) Pivot to wide matrices (gene x cohort) for betas and SEs
# ------------------------
beta_wide <- long %>%
  dplyr::select(gene_id, cohort, log2FoldChange) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = log2FoldChange)

se_wide <- long %>%
  dplyr::select(gene_id, cohort, lfcSE) %>%
  tidyr::pivot_wider(names_from = cohort, values_from = lfcSE)

# Convert to data.table for speed
beta_dt <- as.data.table(beta_wide)
se_dt <- as.data.table(se_wide)

# ------------------------
# 4) Inverse-variance fixed-effect meta-analysis per gene (vectorized)
#    For each gene: weights w_i = 1 / se_i^2
#    beta_meta = sum(w * beta) / sum(w)
#    se_meta = sqrt(1 / sum(w))
#    z_meta = beta_meta / se_meta
#    p_meta = 2 * pnorm(-abs(z_meta))
# ------------------------
cat(crayon::green("Computing inverse-variance fixed-effect meta-analysis per gene...\n"))
meta_fixed <- function(beta_row, se_row) {
  # Input are numeric vectors with NAs for missing cohorts
  valid <- !is.na(beta_row) & !is.na(se_row) & se_row > 0
  if(sum(valid) == 0) return(list(beta_meta=NA, se_meta=NA, z=NA, p=NA, k=0))
  b <- beta_row[valid]
  s <- se_row[valid]
  w <- 1/(s^2)
  beta_meta <- sum(w * b) / sum(w)
  se_meta <- sqrt(1 / sum(w))
  z <- beta_meta / se_meta
  p <- 2 * pnorm(-abs(z))
  return(list(beta_meta=beta_meta, se_meta=se_meta, z=z, p=p, k=sum(valid)))
}

# apply row-wise
n_genes <- nrow(beta_dt)
res_fixed <- vector("list", n_genes)
for(i in seq_len(n_genes)) {
  beta_row <- as.numeric(beta_dt[i, -1, with=FALSE]) # skip gene column
  se_row <- as.numeric(se_dt[i, -1, with=FALSE])
  res_fixed[[i]] <- meta_fixed(beta_row, se_row)
}
res_fixed_df <- do.call(rbind, lapply(res_fixed, function(x) as.data.frame(x)))
res_fixed_df$gene <- beta_dt$gene

# ------------------------
# 5) Random-effects meta-analysis (DerSimonian-Laird) per gene using 'metafor'
#    This step is slower because many small meta-analyses are run (one per gene).
#    We'll run for genes that have >=2 cohorts with non-missing SE.
#    Optionally parallelize using BiocParallel or foreach.
# ------------------------
cat(crayon::green("Computing random-effects meta-analysis (DerSimonian-Laird) per gene using the package 'metafor'...\n"))
# Helper to run metafor::rma for one gene
meta_random_rma <- function(beta_row, se_row) {
  valid <- !is.na(beta_row) & !is.na(se_row) & se_row > 0
  if(sum(valid) < 2) return(list(beta_re=NA, se_re=NA, z_re=NA, p_re=NA, tau2=NA, k=sum(valid)))
  yi <- beta_row[valid]
  sei <- se_row[valid]
  # try-catch: some genes may fail to converge
  fit <- tryCatch({
    rma(yi=yi, sei=sei, method="DL") # DerSimonian-Laird
  }, error = function(e) NULL, warning = function(w) NULL)
  if(is.null(fit)) return(list(beta_re=NA, se_re=NA, z_re=NA, p_re=NA, tau2=NA, k=sum(valid)))
  beta_re <- fit$b[1]
  se_re <- fit$se
  z_re <- beta_re / se_re
  p_re <- 2 * pnorm(-abs(z_re))
  tau2 <- fit$tau2
  return(list(beta_re=beta_re, se_re=se_re, z_re=z_re, p_re=p_re, tau2=tau2, k=sum(valid)))
}

# Option: parallelize
BPPARAM <- MulticoreParam(workers = max(1, parallel::detectCores()-1)) # adjust if on Windows use SnowParam
res_random <- bplapply(seq_len(n_genes), function(i) {
  beta_row <- as.numeric(beta_dt[i, -1, with=FALSE])
  se_row <- as.numeric(se_dt[i, -1, with=FALSE])
  meta_random_rma(beta_row, se_row)
}, BPPARAM = BPPARAM)

res_random_df <- do.call(rbind, lapply(res_random, function(x) as.data.frame(x)))
res_random_df$gene <- beta_dt$gene

# ------------------------
# 6) Combine fixed and random results into one table
# ------------------------
meta_results <- left_join(res_fixed_df, res_random_df, by="gene")
# reorder columns
meta_results <- meta_results %>% dplyr::select(gene, k.x, beta_meta, se_meta, z, p,
                                        beta_re, se_re, z_re, p_re, tau2)
colnames(meta_results)[2] <- "k" # number of cohorts used

# Save meta gene-level table
# fwrite(meta_results, "meta_gene_level_results.csv")
fwrite(meta_results, here(outpath, "meta_gene_level_results.csv"))
cat(crayon::green("Meta-Analysis Done! Moving on to GSEA...\n"))

# ------------------------
# 7) Create ranked gene list for GSEA using meta results
#    Approach: use signed statistic, e.g. signed -log10(p) * sign(beta_meta)
#    Many people use rank by beta_meta, or by z, or by sign(beta)*(-log10(p)).
# ------------------------
# We'll use fixed-effect z if available, else random-effect z, else beta.
meta_results <- meta_results %>%
  mutate(stat_for_ranking = ifelse(!is.na(z), z,
                                   ifelse(!is.na(z_re), z_re, beta_meta))) %>%
  # ensure no NA
  mutate(stat_for_ranking = ifelse(is.na(stat_for_ranking), 0, stat_for_ranking))

ranked <- meta_results %>%
  arrange(desc(stat_for_ranking)) %>%
  dplyr::select(gene, stat_for_ranking)

# Build numeric named vector for fgsea; names = genes
ranks_vec <- ranked$stat_for_ranking
names(ranks_vec) <- ranked$gene

# Offset ranks by just a tiny factor
offs <- (as.integer(factor(names(ranks_vec))) %% 1000) * 1e-12
ranks_vec <- ranks_vec + offs

# ------------------------
# 8) Load gene sets (GMT or msigdbr)
# ------------------------
# Example: use MSigDB Hallmark gene sets via msigdbr (human)
# BiocManager::install("msigdbr")
# library(msigdbr)
# hallmarks <- msigdbr(species = "Homo sapiens", category = "H") %>%
#   split(x = .$gene_symbol, f = .$gs_name)
# hallmarks
# Or load a local GMT:
gmt <- fgsea::gmtPathways(gmt.path)
# gmt <- fgsea::gmtPathways("~/bioinf-tools/db/molsigs/msigdb/h.all.v2024.1.Hs.symbols.gmt")
# gmt <- fgsea::gmtPathways("~/bioinf-tools/db/molsigs/TARSA/TARSA_IFA/TARSA_IFA.gmt")

# For now assume you have a list of pathways called 'pathways_list'

# Example: placeholder (replace with your pathway list)
# pathways_list <- gmt

# ------------------------
# 9) Run fgsea on the ranked vector
# ------------------------

# Use fgseaMultilevel for accurate p-values
# adjust parameters depending on gene set size & desired permutations
fgsea_res <- fgsea::fgseaMultilevel(pathways = gmt,
                                    stats = ranks_vec,
                                    minSize = 15,
                                    maxSize = 2000,
                                    eps = 0)  # eps=0 gives exact p-values

# Tidy results
fgseaResTidy <- fgsea_res %>%
  as.data.frame() %>%
  dplyr::filter(padj < 0.05) %>%
  arrange(padj) # sort by adjusted p-value

# fgseaResTidy <- dplyr::filter(fgseaResTidy, fgseaResTidy$padj < 0.1)

fgseaResTidy$Direction <- recode(factor(fgseaResTidy$NES > 0), "TRUE" = "Upregulated", "FALSE" = "Downregulated")

# Convert leading edge list into vectors for convenience
fgseaResTidy <- fgseaResTidy %>%
  mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ", "))) %>%
  mutate(leadingEdge = as.character(leadingEdge)) # Force character

custom_colors <- c("Upregulated" = "#077f97", "Downregulated" = "#400257")

# title <- design$test
cat(crayon::green("Serialising results and plots...\n"))

p <- ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Direction)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score", title = contrast) + # title=title
  theme(text = element_text(size = 8),
        title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = custom_colors) +
  theme_prism( base_size = 12)

jpeg(here(outpath, 'meta_NES_barplot.png'),
     width = 3200, height = 1800, res = 300)
print(p)
dev.off()

pdf(here(outpath, "meta_NES_barplot.pdf"), width = 12, height = 10)   # large canvas
print(p)
dev.off()

# fwrite(fgsea_res_tidy, "meta_fgsea_results.csv")
fwrite(fgsea_res_tidy, here(outpath, "meta_fgsea_results.csv"))

# ------------------------
# 10) Diagnostics / plots
# ------------------------
# Volcano-like plot for meta-analysis
# library(ggplot2)
# meta_results <- meta_results %>% mutate(logp= -log10(p))
# ggplot(meta_results, aes(x=beta_meta, y=logp)) +
#   geom_point(alpha=0.3) +
#   theme_minimal() +
#   ggtitle("Meta (fixed-effect) volcano: combined log2FC vs -log10(p)")
#
# # Plot top leading edge genes for top pathways
# topPath <- fgsea_res_tidy$pathway[1]
# lead <- fgsea_res_tidy[fgsea_res_tidy$pathway == topPath, "leadingEdge"][[1]]
# head(lead)

# ------------------------
# 11) Notes & optional extras
# ------------------------
# - If your cohorts use different gene identifiers, map them to a common ID before meta-analysis
#   (use biomaRt or org.Hs.eg.db). For mouse-human mixes, map orthologs then aggregate duplicates (e.g. take mean).
# - If you only have GSEA outputs (NES/padj) and NOT gene-level stats, you cannot apply this meta-gene approach;
#   you would need alternative approaches (vote counting, combining pathway p-values via Fisher, etc).
# - You can swap the meta-analysis step for more advanced methods (e.g. Stouffer's Z, or adaptive tests).
# - For very large gene sets / many cohorts, consider more efficient vectorized implementations or C-level code.
# - Save intermediate results and examine genes with high tau2 (heterogeneity) to identify discordant cohorts.
