/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use num_traits::Float;

use crate::{test_report::Case, StrErr};

pub fn parse_wolfram_str<T: Float>(s: &str) -> Result<T, StrErr> {
    Ok(match s.trim() {
        "Infinity" => T::infinity(),
        "-Infinity" => T::neg_infinity(),
        "ComplexInfinity" => T::infinity(),
        "Indeterminate" => T::nan(),
        "NaN" => T::nan(),
        _ => T::from(s.parse::<f64>().map_err(|_| "Cannot parse float")?)
            .expect("Cannot parse float"),
    })
}

/// Reads wolfram data from file and returns a vector of Cases
pub fn read_wolfram_data<T: Float>(file_path: &str) -> Result<Vec<Case<T>>, StrErr> {
    use csv::ReaderBuilder;
    use std::fs::File;
    use std::path::Path;

    let path = Path::new(file_path);
    let file = File::open(path).map_err(|_| {
         eprintln!("Test data not found: {file_path} Download from: https://github.com/p-sira/ellip/tree/main/tests/data");
         "Test data not found."
     })?;
    let mut reader = ReaderBuilder::new().has_headers(false).from_reader(file);

    let mut results = Vec::new();

    for record in reader.records() {
        let row = record.expect("Cannot read row");
        let len = row.len();

        let inputs = row
            .iter()
            .take(len - 1)
            .map(|s| parse_wolfram_str(s).expect(&format!("Cannot parse data {s}")))
            .collect();
        let expected = {
            let s = &row[len - 1];
            parse_wolfram_str(s).expect(&format!("Cannot parse data {s}"))
        };
        results.push(Case { inputs, expected });
    }

    Ok(results)
}
