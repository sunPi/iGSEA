# iGSEA

iGSEA integrates independent differential-expression experiments at the gene
level and then performs one preranked GSEA analysis. It does **not** combine NES
values from separate enrichment analyses.

For every gene, iGSEA pools cohort-level `log2FoldChange` estimates and their
standard errors. It reports both fixed- and random-effects meta-analysis, then
ranks genes using the model selected by `--model`.

## Input

`--data` must be a directory containing at least two CSV files, one per
independent cohort or experiment. Every file requires:

- `gene_id`
- `log2FoldChange`
- `lfcSE`, or `stat` so that `lfcSE = abs(log2FoldChange / stat)` can be derived

All experiments must use the same gene identifier system and the same contrast
direction. For example, every log2 fold change must represent `mutant - wild
type`, never a mixture of that contrast and its reverse.

iGSEA does not deduplicate genes. A duplicated `gene_id` causes a clear error;
run `handleGeneDups` upstream. The former implementation that retained the
first duplicate is preserved at `legacy/main_with_deduplication.R`.

## Models

The default is random effects with restricted maximum likelihood:

```text
--model random --method REML
```

Random effects allows the true gene effect to vary between cohorts and adds the
estimated between-cohort variance, `tau2`, to the inverse-variance weights.
REML is a widely used baseline estimator of `tau2`.

Fixed effects can be selected with:

```text
--model fixed
```

Fixed effects assumes all cohorts estimate one common true effect and that
observed differences arise from cohort-level estimation error alone.

Both models are always calculated and written to the gene-level output. The
selected model only determines the signed Z-statistic supplied to fgsea.

Other supported random-effects estimators are `ML`, `DL`, `HE`, `HS`, `SJ`,
`EB`, and `PM`. `--method DL` reproduces the DerSimonian-Laird choice in the
original script.

## Usage

```bash
iGSEA \
  --data dge_results \
  --gene-set hallmark.gmt \
  --contrast "Mutated vs wild type"
```

Fixed-effect ranking:

```bash
iGSEA -d dge_results -g hallmark.gmt -c "Mutated vs wild type" \
  --model fixed
```

The output directory defaults to `results`. Run `iGSEA --help` for all options.

## Outputs

The existing output names and plot presentation are retained:

- `meta_gene_level_results.csv`
- `meta_fgsea_results.csv` — pathways passing the adjusted p-value threshold
- `meta_fgsea_results_all.csv` — complete enrichment results
- `meta_NES_barplot.png`
- `meta_NES_barplot.pdf`

The gene-level table contains:

- `k`: number of contributing cohorts
- `beta_fixed`, `se_fixed`, `z_fixed`, `p_fixed`, `padj_fixed`
- `beta_random`, `se_random`, `z_random`, `p_random`, `padj_random`
- `tau2`: estimated between-cohort variance
- `Q` and `Q_pvalue`: Cochran's heterogeneity test
- `I2`: estimated percentage of observed variation attributable to
  heterogeneity rather than sampling error
- `direction_agreement`: fraction of contributing cohort effects agreeing with
  the majority direction

Cochran's Q and I2 describe heterogeneity; they are not used as the standard
directional GSEA ranking. The selected fixed- or random-effects Z-statistic is
used because it retains both effect direction and strength of evidence.

Positive pathways use teal (`#077F97`) and negative pathways use purple
(`#400257`).

## Dependencies

Required R packages:

- `data.table`
- `dplyr`
- `metafor`
- `fgsea`
- `BiocParallel`
- `ggplot2`
- `ggprism`

iGSEA no longer requires the private `pipecore` package.

## Limitations

- Random-effects heterogeneity is difficult to estimate with very few cohorts.
- Meta-analysis cannot determine whether heterogeneity is caused by biology,
  batch, platform, treatment history, or another study-level difference.
- Cohorts are assumed to be independent.
- Effect estimates must be scientifically comparable across experiments.
