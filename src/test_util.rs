/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

pub const RTOL: f64 = 5.0 * f64::EPSILON;

#[macro_export]
macro_rules! compare_test_data {
    ($file_path: expr, $func: expr, $rtol: expr) => {
        use std::fs::File;
        use std::io::{self, BufRead};
        use std::path::Path;

        // Check if the file exists
        let path = Path::new($file_path);
        if !path.exists() {
            panic!("Test data not found: {}", $file_path);
        }

        let file = File::open($file_path).expect("Cannot open file");
        let reader = io::BufReader::new(file);

        for (line_number, line) in reader.lines().enumerate() {
            let line = line.expect("Cannot read line");

            let parts: Vec<&str> = line.split_whitespace().collect();
            let inputs: Vec<f64> = parts[..parts.len() - 1]
                .iter()
                .map(|&v| v.parse().expect("Cannot parse input(s) as a number"))
                .collect();

            let expected: f64 = parts[parts.len() - 1]
                .parse()
                .expect("Cannot parse expected value as a number");

            let result = $func(&inputs);
            // Formula for rtol
            if (result - expected).abs() > $rtol * expected {
                panic!(
                    "Test failed on line {}: input = {:?}, expected = {}, got = {}",
                    line_number + 1,
                    inputs,
                    expected,
                    result
                );
            }
        }
    };
}
