use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use serde::Deserialize;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

pub const N_PCS: usize = 30;

/// One SNP entry from glad_pca_weights.tsv.gz
#[derive(Debug, Clone)]
pub struct SnpWeight {
    pub chrom: String,
    pub pos: u64,
    pub effect_allele: String,
    pub other_allele: String,
    pub weights: [f64; N_PCS],
}

/// Contents of glad_meta.json
#[derive(Debug, Deserialize)]
pub struct GladMeta {
    pub reference_build: String,
    pub n_pcs: usize,
    pub n_snps: usize,
    pub age_mean: f64,
    pub age_sd: f64,
}

pub fn load_meta(path: &Path) -> Result<GladMeta> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let meta: GladMeta = serde_json::from_reader(BufReader::new(file))
        .with_context(|| format!("parsing {}", path.display()))?;
    Ok(meta)
}

pub fn load_weights(path: &Path) -> Result<Vec<SnpWeight>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let decoder = GzDecoder::new(file);
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(decoder);

    let headers = rdr.headers()?.clone();
    // Expect: chrom, pos, effect_allele, other_allele, pc1..pc30
    let n_pc_cols = headers.len().saturating_sub(4);
    anyhow::ensure!(
        n_pc_cols == N_PCS,
        "expected {} PC columns in weights file, found {}",
        N_PCS,
        n_pc_cols
    );

    let mut snps = Vec::new();
    for (i, result) in rdr.records().enumerate() {
        let rec = result.with_context(|| format!("reading weights row {}", i + 1))?;
        let chrom = rec[0].to_string();
        let pos: u64 = rec[1]
            .parse()
            .with_context(|| format!("parsing pos at row {}", i + 1))?;
        let effect_allele = rec[2].to_string();
        let other_allele = rec[3].to_string();
        let mut weights = [0f64; N_PCS];
        for j in 0..N_PCS {
            weights[j] = rec[4 + j]
                .parse()
                .with_context(|| format!("parsing pc{} at row {}", j + 1, i + 1))?;
        }
        snps.push(SnpWeight { chrom, pos, effect_allele, other_allele, weights });
    }
    Ok(snps)
}
