# glad-prep

Prepare genomic queries for control matching at [glad.igs.umaryland.edu](https://glad.igs.umaryland.edu).

## Overview

glad-prep is a command-line tool that prepares your cohort's genomic data for submission to the GLAD control matching service. Given an imputed VCF, it extracts allele counts at a set of reference SNPs and projects your samples into a shared PC space, then summarizes the resulting ancestry distribution as a Gaussian mixture model. The output is a small archive that you upload to the GLAD website, where it is used to find well-matched controls from the GLAD database.

No individual-level genotype data leaves your machine. The archive contains only summary statistics (allele counts) and a fitted distribution over PC space — neither of which can be used to recover individual genotypes.

## Prerequisites

- A genome-wide imputed VCF (GRCh38, bgzipped) with an accompanying `.tbi` tabix index. We recommend the [University of Michigan Imputation Server](https://imputationserver.sph.umich.edu) if you need to impute your data.
- Reference files downloaded from [glad.igs.umaryland.edu](https://glad.igs.umaryland.edu): `glad_pca_weights.tsv.gz` and `glad_meta.json`.
- *(Optional)* A sample metadata TSV with sex and/or age for your samples (see [Sample metadata](#sample-metadata)). Including metadata enables sex-stratified and age-aware matching.

## Installation

Pre-built binaries for Linux and macOS are available on the [GLAD website](https://glad.igs.umaryland.edu).

To compile from source, install [Rust](https://rustup.rs) and run:

```sh
cargo build --release
# binary is at target/release/glad-prep
```

## Quick start

First, run `check` to confirm your VCF has adequate coverage of the reference SNPs:

```sh
glad-prep check --vcf your_cohort.vcf.gz
```

Then run `prepare` to generate the submission archive:

```sh
glad-prep prepare --vcf your_cohort.vcf.gz --sample-meta metadata.tsv
```

This writes `query.glad.gz` in the current directory. Upload it at [glad.igs.umaryland.edu](https://glad.igs.umaryland.edu).

Both commands expect `glad_pca_weights.tsv.gz` and `glad_meta.json` in the current directory (the defaults). If you place them elsewhere, use `--weights` and `--meta` to point to them explicitly.

## Command reference

### `check`

Validates that your VCF contains the reference SNPs and reports coverage. Run this before `prepare` to catch problems early.

```
glad-prep check --vcf <path> [--weights <path>]
```

### `prepare`

Extracts allele counts, projects samples to PC space, fits a Gaussian mixture model, and writes the submission archive.

```
glad-prep prepare --vcf <path> [--weights <path>] [--meta <path>]
                  [--sample-meta <path>] [--output <path>]
```

| Flag | Default | Description |
|---|---|---|
| `--vcf` | *(required)* | Imputed VCF.gz with .tbi index |
| `--weights` | `glad_pca_weights.tsv.gz` | GLAD PCA weights file |
| `--meta` | `glad_meta.json` | GLAD reference metadata |
| `--sample-meta` | *(none)* | Sample metadata TSV (sex, age) |
| `--output` | `query.glad.gz` | Output archive path |

## Sample metadata

The `--sample-meta` file is a tab-separated file with the following columns:

```
sample_id	sex	age
SAMPLE_001	1	54.2
SAMPLE_002	0	61.8
```

- `sample_id` must match the sample IDs in the VCF header exactly.
- `sex`: `1` = male, `0` = female.
- `age`: continuous value in years.

All three columns are optional, but if `sex` or `age` is present it must be present for every sample — partial data is not supported. If `sex` is absent or incomplete, a single ancestry distribution is fitted across all samples. If `age` is absent or incomplete, age is excluded from the distribution. Providing both enables the most precise matching.

The minimum supported cohort size is 100 samples.

## Output

`query.glad.gz` is a gzip-compressed JSON file containing:

- **Allele counts** at each reference SNP (alt allele count and total observed alleles), used for genomic inflation (lambda) calculation during matching.
- **Ancestry distribution** — a Gaussian mixture model fitted in PC space, used by the sinkhorn-based matching algorithm to find controls with a similar ancestry profile.

Upload this file at [glad.igs.umaryland.edu](https://glad.igs.umaryland.edu).

## Privacy

glad-prep is designed so that individual-level data never needs to leave your institution. The allele counts in the output are population-level summary statistics. The ancestry distribution is a parametric summary (mixture model means, covariances, and weights) rather than a record of individual PC coordinates. Both are safe to share under standard data sharing practices.
