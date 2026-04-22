use crate::reference::{SnpWeight, N_PCS};
use crate::vcf::VcfData;
use ndarray::Array2;

/// Project samples into PC space.
/// Returns an (n_samples × N_PCS) matrix.
/// SNPs absent from the VCF contribute 0 to all projections.
pub fn project(vcf_data: &VcfData, snps: &[SnpWeight]) -> Array2<f64> {
    let n_samples = vcf_data.samples.len();
    let mut pc_coords = Array2::<f64>::zeros((n_samples, N_PCS));

    for (snp_idx, snp_data_opt) in vcf_data.snps.iter().enumerate() {
        let snp_data = match snp_data_opt {
            Some(d) => d,
            None => continue,
        };
        let weights = &snps[snp_idx].weights;
        for (sample_idx, &dosage) in snp_data.dosages.iter().enumerate() {
            if dosage == 0 {
                continue;
            }
            let d = dosage as f64;
            for pc in 0..N_PCS {
                pc_coords[(sample_idx, pc)] += d * weights[pc];
            }
        }
    }

    pc_coords
}
