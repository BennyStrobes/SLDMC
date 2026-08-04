# SLDMC

SLDMC estimates the correlation between *predicted* variant effects on gene expression and *true* causal cis-eQTL
effect sizes, stratified by variant-gene annotations. A predicted variant-gene effect is, for example,
the effects predicted by a sequence-to-expression model such as Borzoi ir AlphaGenome.

The estimator requires only marginal eQTL summary statistics and LD; it does not require individual-level
expression data, and it does not require knowing which variants are causal. Because the correlation is
estimated separately within each category of an annotation, SLDMC can be used to ask where a predictor
is well correlated with the truth and where it is not (e.g. by distance to TSS, by allele frequency, by predicted effect
magnitude, or within functional annotations). It can also be used to assess calibration.

## Method overview

Analysis is performed across genes, and across all the cis variants of that gene. For each category `c` of an annotation, SLDMC estimates three quantities:

1. **Calibration slope** `α_c` — the slope from regressing marginal eQTL effect sizes on LD-weighted
   predicted effects, restricted to variants in category `c`.
2. **Predicted-effect variance** `σ²_c` — the average squared predicted effect of variants in category `c`.
3. **Per-SNP cis-eQTL heritability** `τ_c` — obtained by regressing `β̂² − se²` on stratified LD scores
   (a stratified LD score regression, using the unbiased squared-LD estimator
   `r² − (1 − r²)/(N − 2)`).

The correlation for the category is then

```
correlation_c = α_c * sqrt(σ²_c / τ_c)
```

Both eQTL effect sizes and predicted effect sizes are put on a standardized-genotype scale before the
regressions are run. LD is computed in-sample, from the genotypes of the individuals in the eQTL study.

Standard errors and confidence intervals come from a bootstrap that resamples **genes** with replacement
(100 iterations). Each gene contributes a small set of sufficient statistics (`XᵀX` and `Xᵀy` per
regression), so bootstrap iterations only require re-summing per-gene matrices rather than re-reading
genotypes.

## Installation

```bash
git clone https://github.com/BennyStrobes/SLDMC
cd SLDMC
conda env create -f environment.yml
conda activate sldmc
```

The environment pins `python=3.12`, `numpy=2.3.4`, and `pandas-plink=2.3.2`. Python 3.12 is the newest
interpreter for which `pandas-plink` 2.3.2 packages are available on conda-forge.

## Running SLDMC

```bash
python sldmc.py \
    --est-borzoi-effect-size-file predicted_effects.txt.gz \
    --est-eqtl-effect-size-file eqtl_sumstats.txt.gz \
    --sim-variant-gene-annotation-file annotations.txt.gz \
    --genotype-plink-filestem /path/to/genotype_chr \
    --genotype-sample-mapping-file genotype_sample_mapping.txt \
    --ld-corr-output-stem /path/to/output/my_analysis \
    --bootstrapped-gene-set-filestem /path/to/gene_sets/bootstrap_ \
    --n-bootstraps 500 \
    --weighted True
```

| Argument | Required | Description |
| --- | --- | --- |
| `--est-borzoi-effect-size-file` | yes | Predicted variant-gene effect sizes |
| `--est-eqtl-effect-size-file` | yes | Marginal cis-eQTL summary statistics |
| `--sim-variant-gene-annotation-file` | yes | Annotation category of each variant-gene pair |
| `--genotype-plink-filestem` | yes | PLINK filestem, chromosome number appended |
| `--genotype-sample-mapping-file` | yes | Which genotype samples make up the eQTL sample |
| `--ld-corr-output-stem` | yes | Output filestem |
| `--bootstrapped-gene-set-filestem` | no | Pre-specified bootstrap gene sets (see below) |
| `--n-bootstraps` | no | Number of bootstrap iterations. Default `500` |
| `--weighted` | no | Use LD score regression weights. Default `True` |

A gene is analyzed only if it appears in all three of the effect, eQTL, and annotation files, and has at
least 10 cis variants that survive the allele checks described below. All 22 autosomes are processed.

## Input file formats

All three of the variant-gene files are gzipped, tab-delimited, and carry a header line. Their first six
columns are the same: gene ID, variant ID, chromosome, position, allele 0, allele 1. Variant IDs must be
consistent across files and with the PLINK `.bim`.

### Predicted effect size file (`--est-borzoi-effect-size-file`)

| Column | Name | Description |
| --- | --- | --- |
| 1 | `gene` | Gene ID |
| 2 | `variant` | Variant ID |
| 3 | `chr` | Chromosome (1-22, no `chr` prefix) |
| 4 | `snp_pos` | Position |
| 5 | `a0` | Allele 0 |
| 6 | `a1` | Allele 1 |
| 7 | `borzoi_effect_size` | Predicted effect of `a1` relative to `a0` |

Predicted effects are supplied on the raw genotype scale; SLDMC multiplies them by the genotype standard
deviation internally. Both the standardized and unstandardized variances are reported in the output.

### eQTL summary statistics file (`--est-eqtl-effect-size-file`)

| Column | Name | Description |
| --- | --- | --- |
| 1 | `gene` | Gene ID |
| 2 | `variant` | Variant ID |
| 3 | `chr` | Chromosome |
| 4 | `snp_pos` | Position |
| 5 | `A0` | Allele 0 |
| 6 | `A1` | Allele 1 |
| 7 | `eqtl_effect_size` | Marginal eQTL effect size (slope) |
| 8 | `eqtl_effect_size_se` | Standard error of the effect size |
| 9 | `N` | eQTL sample size |
| 10 | `maf` | Minor allele frequency |
| 11 | `genotype_sdev` | *(optional)* Genotype standard deviation |

If an 11th column is present it is used as the per-variant genotype standard deviation. If it is absent,
or its value is the literal string `nan`, the standard deviation is computed from the MAF as
`sqrt(2 * maf * (1 - maf))`. This lets a file mix variants that have a measured standard deviation with
variants that do not.

### Variant-gene annotation file (`--sim-variant-gene-annotation-file`)

After the six shared columns, this file has **one column per annotation**. Each entry is the 0-based
index of the category that variant-gene pair falls into for that annotation, or `-1` if the pair falls
into no category of that annotation. The column headers are the annotation names.

```
gene             variant               chr  snp_pos  a0  a1  intercept  dist_to_tss_bins  af_bins
ENSG00000187634  chr1_925952_G_A_b38   1    925952   G   A   0          2                 1
ENSG00000187634  chr1_926431_C_T_b38   1    926431   C   T   0          3                 -1
```

### Annotation category file

SLDMC also reads a companion file describing the categories of each annotation. Its path is derived from
the annotation file by replacing the `.txt.gz` suffix with `_categories.txt`; it is plain text,
tab-delimited, with a header.

| Column | Name | Description |
| --- | --- | --- |
| 1 | `anno_name` | Annotation name, matching a column of the annotation file |
| 2 | `source` | Free-text label describing where the annotation came from |
| 3 | `category_index` | 0-based category index |
| 4 | `category_name` | Name of that category |

Every annotation column in the annotation file must appear here, and each annotation's categories must
be listed in increasing `category_index` order starting from 0.

### Genotype files (`--genotype-plink-filestem`)

Per-chromosome PLINK `.bed`/`.bim`/`.fam` files, given as a filestem with the chromosome number appended.
For example `--genotype-plink-filestem /path/to/genotype_chr` expects `/path/to/genotype_chr1.bed`
through `/path/to/genotype_chr22.bed`.

Variants are matched to the other input files by the `.bim` variant ID. A variant is dropped if its
`.bim` alleles do not match, as an unordered pair, the alleles listed in the effect and summary statistic
files. Effect signs in the other files are flipped where needed to match the `.bim` allele orientation,
so the input files do not need to agree on which allele is the reference. Missing genotypes are
mean-imputed before LD is computed.

### Genotype sample mapping file (`--genotype-sample-mapping-file`)

Plain text, no header, one integer per line: the 0-based index of a sample in the PLINK `.fam` ordering.
Together these lines identify the individuals in the eQTL study, so that LD is computed in-sample rather
than from an external reference panel.

### Bootstrapped gene sets (`--bootstrapped-gene-set-filestem`, optional)

If omitted, each bootstrap iteration resamples genes independently within the run. If supplied, SLDMC
reads `<filestem>0.txt` through `<filestem><n-1>.txt`, where `n` is `--n-bootstraps`; each file has a
header line followed by one gene ID per line, sampled with replacement. Gene IDs not present in the
current analysis are dropped.

Supplying pre-specified gene sets is what makes results comparable across separate SLDMC runs: when
every run resamples the same genes on iteration `i`, the bootstrap distributions are aligned and can be
combined across runs, for example when meta-analyzing across tissues.

## Output files

Two tab-delimited files are written, both prefixed with `--ld-corr-output-stem`.

### `<stem>_bootstrap_summary.txt`

One row per (annotation, category, sample). The `observed` row (with `sample_iter` of `-1`) holds the
point estimate; the remaining rows hold the `--n-bootstraps` bootstrap iterations.

| Column | Description |
| --- | --- |
| `annotation_name` | Annotation name |
| `category_name` | Category within that annotation |
| `sample_label` | `observed` or `bootstrap` |
| `sample_iter` | Bootstrap iteration, or `-1` for the observed estimate |
| `calibration_slope` | Calibration slope `α_c` |
| `per_snp_eqtl_h2` | Per-SNP cis-eQTL heritability `τ_c` |
| `borzoi_variance` | Variance of standardized predicted effects `σ²_c` |
| `unstandardized_borzoi_variance` | Variance of raw predicted effects |
| `correlation` | `α_c * sqrt(σ²_c / τ_c)` |

`correlation` is `nan` where `σ²_c / τ_c` is negative or non-finite, which can happen for categories
with a small or noisily estimated heritability.

### `<stem>_bootstrap_stats.txt`

The bootstrap distribution of each quantity above, summarized into one row per
(annotation, category, quantity).

| Column | Description |
| --- | --- |
| `annotation_name` | Annotation name |
| `category_name` | Category within that annotation |
| `output_name` | Which quantity this row summarizes |
| `mean` | The observed (non-bootstrap) point estimate |
| `bootstrapped_mean` | Mean across bootstrap iterations |
| `bootstrap_se` | Standard deviation across bootstrap iterations |
| `gaussian_z_score` | `mean / bootstrap_se` |
| `empirical_ci_lower` | 2.5th percentile of the bootstrap distribution |
| `empirical_ci_upper` | 97.5th percentile of the bootstrap distribution |

## Compute requirements

SLDMC holds one chromosome of genotypes in memory at a time and computes an LD matrix per gene. In the
analyses reported in the paper, each run over a genome-wide set of genes was submitted with a 53GB memory
allocation and a four hour wall-clock limit.

## Citation

*Manuscript in preparation.*

## Contact

Ben Strober — bstrober3@gmail.com. Questions and bug reports are welcome via
[GitHub issues](https://github.com/BennyStrobes/SLDMC/issues).
