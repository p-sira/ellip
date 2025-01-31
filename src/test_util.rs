/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

pub const RTOL: f64 = 5.0 * f64::EPSILON;

#[macro_export]
macro_rules! compare_test_data {
    ($file_path:expr, $func:expr, $rtol:expr) => {
        compare_test_data!($file_path, $func, $rtol, 0.0)
    };
    ($file_path:expr, $func:expr, $rtol:expr, $atol:expr) => {
        {
            use std::fs::File;
            use std::io::{self, BufRead};
            use std::path::Path;

            let path = Path::new($file_path);
            if !path.exists() {
                panic!("Test data not found: {}", $file_path);
            }

            let file = File::open($file_path).expect("Cannot open file");
            let reader = io::BufReader::new(file);

            let mut test_fail = false;
            let mut worst_line = 0;
            let mut worst_params: Vec<f64> = Vec::new();
            let mut worst_errors = [0.0; 5];
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
                let error = (result - expected).abs();
                let tol = $atol + $rtol * expected.abs();
                if error > tol {
                    let rel = error / expected.abs();
                    test_fail = true;
                    if rel > worst_errors[1] {
                        worst_line = line_number + 1;
                        worst_errors = [error, rel, tol, expected, result];
                        worst_params = inputs.clone();
                    }
                    eprintln!(
                        "Test failed on line {}: input = {:?}, expected = {:?}, got = {:?} \n error = {:?}, rel = {:?}, abs = {:?}, rtol = {:?}, atol = {:?}\n",
                        line_number + 1,
                        inputs,
                        expected,
                        result,
                        error,
                        rel,
                        tol,
                        $rtol,
                        $atol,
                    );
                }
            }

            if test_fail {
                panic!(
                    "Worst on line {}: input = {:?}, expected = {:?}, got = {:?} \n error = {:?}, rel = {:?}, abs = {:?}, rtol = {:?}, atol = {:?}",
                    worst_line,
                    worst_params,
                    worst_errors[3],
                    worst_errors[4],
                    worst_errors[0],
                    worst_errors[1],
                    worst_errors[2],
                    $rtol,
                    $atol,
                );
            }
        }
    };
}

#[macro_export]
macro_rules! assert_close {
    ($expected: expr, $actual: expr, $rtol: expr) => {
        let relative = ($actual - $expected).abs() / $expected;
        if relative > $rtol {
            panic!(
                "Assertion failed: expected = {:?}, got = {:?}, relative = {:?}, rtol = {:?}",
                $expected, $actual, relative, $rtol
            )
        }
    };
}

pub fn linspace(start: f64, end: f64, num: usize) -> Vec<f64> {
    if num < 2 {
        return vec![start];
    }

    let step = (end - start) / (num - 1) as f64;
    let mut result = Vec::with_capacity(num);

    for i in 0..num {
        result.push(start + step * i as f64);
    }

    result
}
