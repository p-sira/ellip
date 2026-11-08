/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use std::path::{Path, PathBuf};
use std::sync::LazyLock;

use crate::StrErr;

/// Resolves Criterion's output directory once, using the same precedence as Criterion.
pub fn criterion_output_directory() -> &'static Path {
    static DIRECTORY: LazyLock<PathBuf> = LazyLock::new(|| {
        if let Some(path) = std::env::var_os("CRITERION_HOME") {
            return PathBuf::from(path);
        }
        if let Some(path) = std::env::var_os("CARGO_TARGET_DIR") {
            return PathBuf::from(path).join("criterion");
        }

        #[derive(serde::Deserialize)]
        struct Metadata {
            target_directory: PathBuf,
        }

        let target_directory = (|| -> Result<PathBuf, Box<dyn std::error::Error>> {
            let output = std::process::Command::new(
                std::env::var_os("CARGO").unwrap_or_else(|| "cargo".into()),
            )
            .args(["metadata", "--no-deps", "--format-version", "1"])
            .output()?;
            if !output.status.success() {
                return Err(String::from_utf8_lossy(&output.stderr).into_owned().into());
            }
            Ok(serde_json::from_slice::<Metadata>(&output.stdout)?.target_directory)
        })();

        match target_directory {
            Ok(path) => path.join("criterion"),
            Err(error) => {
                eprintln!("Cannot resolve Cargo target directory: {error}; using target/criterion");
                PathBuf::from("target/criterion")
            }
        }
    });
    &DIRECTORY
}

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

        let mean = extract_criterion_mean(&estimates_path).unwrap();
        assert_eq!(mean, 150.5);
    }

    #[test]
    fn test_extract_criterion_means() {
        // Create temporary directories with mock estimates.json files
        let temp_dir1 = tempdir().unwrap();
        let temp_dir2 = tempdir().unwrap();

        let estimates_content1 = r#"{"mean":{"point_estimate":100.0}}"#;
        let estimates_content2 = r#"{"mean":{"point_estimate":200.0}}"#;

        let paths = vec![
            temp_dir1.path().join("estimates.json"),
            temp_dir2.path().join("estimates.json"),
        ];

        fs::write(paths[0].clone(), estimates_content1).unwrap();
        fs::write(paths[1].clone(), estimates_content2).unwrap();

        let means = extract_criterion_means(&paths).unwrap();
        assert_eq!(means, vec![100.0, 200.0]);
    }
}
