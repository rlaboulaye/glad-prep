use crate::gmm::Distributions;
use crate::reference::{GladMeta, SnpWeight};
use crate::vcf::VcfData;
use anyhow::{Context, Result};
use flate2::{write::GzEncoder, Compression};
use serde::Serialize;
use std::path::Path;

#[derive(Serialize)]
struct SnpCount {
    chrom: String,
    pos: u64,
    effect_allele: String,
    other_allele: String,
    alt_count: u32,
    n_alleles: u32,
}

#[derive(Serialize)]
struct Output {
    version: &'static str,
    reference_build: String,
    n_samples: usize,
    n_snps_attempted: usize,
    n_snps_found: usize,
    counts: Vec<SnpCount>,
    distributions: Distributions,
}

pub fn write(
    path: &Path,
    snps: &[SnpWeight],
    vcf_data: &VcfData,
    distributions: Distributions,
    glad_meta: &GladMeta,
) -> Result<()> {
    let n_samples = vcf_data.samples.len();
    let n_snps_found = vcf_data.snps.iter().filter(|s| s.is_some()).count();

    let counts: Vec<SnpCount> = snps
        .iter()
        .zip(vcf_data.snps.iter())
        .filter_map(|(snp, data)| {
            data.as_ref().map(|d| SnpCount {
                chrom: snp.chrom.clone(),
                pos: snp.pos,
                effect_allele: snp.effect_allele.clone(),
                other_allele: snp.other_allele.clone(),
                alt_count: d.alt_count,
                n_alleles: d.n_alleles,
            })
        })
        .collect();

    let out = Output {
        version: "1.0",
        reference_build: glad_meta.reference_build.clone(),
        n_samples,
        n_snps_attempted: snps.len(),
        n_snps_found,
        counts,
        distributions,
    };

    let file = std::fs::File::create(path)
        .with_context(|| format!("creating output file {}", path.display()))?;
    let mut gz = GzEncoder::new(file, Compression::default());
    serde_json::to_writer(&mut gz, &out).context("serializing output")?;
    gz.finish().context("finalizing gzip output")?;

    Ok(())
}
