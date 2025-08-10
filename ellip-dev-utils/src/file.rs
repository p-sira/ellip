/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use std::path::PathBuf;

/// Finds test data files matching the pattern: function_name_*.csv
pub fn find_test_files(function_name: &str, dir: &str) -> Vec<PathBuf> {
    // Also try the wolfram directory specifically as fallback
    let pattern = format!("tests/data/{dir}/{function_name}_*.csv");
    glob::glob(&pattern)
        .unwrap()
        .map(|path| path.unwrap())
        .collect()
}
