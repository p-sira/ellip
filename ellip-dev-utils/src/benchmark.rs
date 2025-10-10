/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use std::path::PathBuf;

use crate::StrErr;

/// Criterion benchmark estimates structure
#[derive(Debug, serde::Deserialize)]
pub struct CriterionEstimates {
    pub mean: CriterionMean,
}

/// Criterion mean structure
#[derive(Debug, serde::Deserialize)]
pub struct CriterionMean {
    pub point_estimate: f64,
}

pub fn extract_criterion_mean(path: &PathBuf) -> Result<f64, StrErr> {
    use std::fs;
    let content = fs::read_to_string(path).map_err(|_| "Cannot read estimates.json file")?;

    let estimates: CriterionEstimates =
        serde_json::from_str(&content).map_err(|_| "Cannot parse estimates.json file")?;

    Ok(estimates.mean.point_estimate)
}

pub fn extract_criterion_means(paths: &[PathBuf]) -> Result<Vec<f64>, StrErr> {
    paths
        .iter()
        .map(|path| extract_criterion_mean(path))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::tempdir;

    #[test]
    fn test_extract_criterion_mean() {
        // Create a temporary directory with a mock estimates.json
        let temp_dir = tempdir().unwrap();
        let estimates_content = r#"{
            "mean": {
                "confidence_interval": {
                    "confidence_level": 0.95,
                    "lower_bound": 100.0,
                    "upper_bound": 200.0
                },
                "point_estimate": 150.5,
                "standard_error": 25.0
            }
        }"#;

        let estimates_path = temp_dir.path().join("estimates.json");
        fs::write(&estimates_path, estimates_content).unwrap();

        let mean = extract_criterion_mean(temp_dir.path().to_str().unwrap()).unwrap();
        assert_eq!(mean, 150.5);
    }

    #[test]
    fn test_extract_criterion_means() {
        // Create temporary directories with mock estimates.json files
        let temp_dir1 = tempdir().unwrap();
        let temp_dir2 = tempdir().unwrap();

        let estimates_content1 = r#"{"mean":{"point_estimate":100.0}}"#;
        let estimates_content2 = r#"{"mean":{"point_estimate":200.0}}"#;

        fs::write(temp_dir1.path().join("estimates.json"), estimates_content1).unwrap();
        fs::write(temp_dir2.path().join("estimates.json"), estimates_content2).unwrap();

        let paths = vec![
            temp_dir1.path().to_str().unwrap(),
            temp_dir2.path().to_str().unwrap(),
        ];

        let means = extract_criterion_means(&paths).unwrap();
        assert_eq!(means, vec![100.0, 200.0]);
    }
}
