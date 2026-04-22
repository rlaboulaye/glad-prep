use crate::reference::GladMeta;
use crate::sample_meta::{MetaMode, SampleMeta, Sex};
use anyhow::{Context, Result};
use linfa::traits::Fit;
use linfa::DatasetBase;
use linfa_clustering::GaussianMixtureModel;
use ndarray::{Array2, Axis};
use serde::Serialize;

const MAX_COMPONENTS: usize = 8;
const GMM_MAX_ITER: u64 = 200;
const GMM_RUNS: usize = 5;

#[derive(Debug, Serialize)]
pub struct FittedGmm {
    pub n_components: usize,
    pub weights: Vec<f64>,
    pub means: Vec<Vec<f64>>,
    /// Full covariance matrices: n_components × n_dims × n_dims
    pub covariances: Vec<Vec<Vec<f64>>>,
}

#[derive(Debug, Serialize)]
pub struct Distributions {
    pub mode: String,
    pub n_dims: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub all: Option<FittedGmm>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub female: Option<FittedGmm>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub male: Option<FittedGmm>,
}

pub fn fit(
    pc_coords: &Array2<f64>,
    sample_order: &[String],
    meta: Option<&SampleMeta>,
    glad_meta: &GladMeta,
) -> Result<Distributions> {
    let mode = meta.map(|m| m.mode).unwrap_or(MetaMode::None);
    let use_sex = matches!(mode, MetaMode::SexAndAge | MetaMode::SexOnly);
    let use_age = matches!(mode, MetaMode::SexAndAge | MetaMode::AgeOnly);

    let mode_str = match mode {
        MetaMode::SexAndAge => "sex_and_age",
        MetaMode::SexOnly => "sex_only",
        MetaMode::AgeOnly => "age_only",
        MetaMode::None => "none",
    }
    .to_string();

    let features = build_features(pc_coords, sample_order, meta, glad_meta, use_age);
    let n_dims = features.ncols();

    if use_sex {
        let meta = meta.unwrap();
        let mut female_rows = Vec::new();
        let mut male_rows = Vec::new();
        for (i, id) in sample_order.iter().enumerate() {
            if meta.map.get(id).map(|s| s.sex == Sex::Female).unwrap_or(true) {
                female_rows.push(i);
            } else {
                male_rows.push(i);
            }
        }

        let female_features = select_rows(&features, &female_rows);
        let male_features = select_rows(&features, &male_rows);

        Ok(Distributions {
            mode: mode_str,
            n_dims,
            all: None,
            female: Some(fit_gmm(&female_features).context("fitting female GMM")?),
            male: Some(fit_gmm(&male_features).context("fitting male GMM")?),
        })
    } else {
        Ok(Distributions {
            mode: mode_str,
            n_dims,
            all: Some(fit_gmm(&features).context("fitting GMM")?),
            female: None,
            male: None,
        })
    }
}

fn build_features(
    pc_coords: &Array2<f64>,
    sample_order: &[String],
    meta: Option<&SampleMeta>,
    glad_meta: &GladMeta,
    use_age: bool,
) -> Array2<f64> {
    let n = sample_order.len();
    let n_pcs = pc_coords.ncols();
    let n_cols = if use_age { n_pcs + 1 } else { n_pcs };
    let mut features = Array2::<f64>::zeros((n, n_cols));

    for i in 0..n {
        for j in 0..n_pcs {
            features[(i, j)] = pc_coords[(i, j)];
        }
    }

    if use_age {
        let meta = meta.unwrap();
        for (i, id) in sample_order.iter().enumerate() {
            if let Some(info) = meta.map.get(id) {
                features[(i, n_pcs)] = (info.age - glad_meta.age_mean) / glad_meta.age_sd;
            }
        }
    }

    features
}

fn select_rows(matrix: &Array2<f64>, rows: &[usize]) -> Array2<f64> {
    let n = rows.len();
    let d = matrix.ncols();
    let mut result = Array2::<f64>::zeros((n, d));
    for (new_i, &old_i) in rows.iter().enumerate() {
        result.row_mut(new_i).assign(&matrix.row(old_i));
    }
    result
}

fn fit_gmm(data: &Array2<f64>) -> Result<FittedGmm> {
    let n = data.nrows();
    let n_dims = data.ncols();
    anyhow::ensure!(n >= 10, "too few samples to fit GMM ({})", n);

    // K heuristic: roughly one component per 200 samples, capped at MAX_COMPONENTS
    let k = (n / 200).clamp(1, MAX_COMPONENTS);

    let dataset = DatasetBase::from(data.clone());

    // Multiple restarts; keep the first successful fit
    // TODO: implement BIC-based K selection with log-likelihood scoring
    let mut last_err = None;
    for _ in 0..GMM_RUNS {
        match GaussianMixtureModel::params(k)
            .max_n_iterations(GMM_MAX_ITER)
            .fit(&dataset)
        {
            Ok(model) => return serialize_gmm(&model, k, n_dims),
            Err(e) => last_err = Some(e),
        }
    }

    Err(anyhow::anyhow!(
        "GMM fitting failed after {} runs: {:?}",
        GMM_RUNS,
        last_err
    ))
}

fn serialize_gmm(
    gmm: &GaussianMixtureModel<f64>,
    k: usize,
    n_dims: usize,
) -> Result<FittedGmm> {
    let weights = gmm.weights().to_vec();
    let means_arr = gmm.means();
    let covs_arr = gmm.covariances();

    let mut means = Vec::with_capacity(k);
    let mut covariances = Vec::with_capacity(k);

    for i in 0..k {
        means.push(means_arr.row(i).to_vec());
        let cov = covs_arr.index_axis(Axis(0), i);
        let cov_rows: Vec<Vec<f64>> = (0..n_dims).map(|j| cov.row(j).to_vec()).collect();
        covariances.push(cov_rows);
    }

    Ok(FittedGmm { n_components: k, weights, means, covariances })
}
