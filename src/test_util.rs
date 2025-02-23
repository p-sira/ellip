/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

#[macro_export]
macro_rules! compare_test_data {
    ($file_path:expr, $func:expr, $rtol:expr) => {
        compare_test_data!($file_path, $func, f64, $rtol, 0.0)
    };
    ($file_path:expr, $func:expr, $t: ident, $rtol:expr) => {
        compare_test_data!($file_path, $func, f64, $rtol, 0.0)
    };
    ($file_path:expr, $func:expr, $rtol:expr, $atol:expr) => {
        compare_test_data!($file_path, $func, f64, $rtol, $atol)
    };
    ($file_path:expr, $func:expr, $t: ident, $rtol:expr, $atol:expr) => {
        {
            use std::fs::File;
            use std::io::{BufRead, BufReader};
            use std::path::Path;
            use std::str::FromStr;

            let path = Path::new($file_path);
            if !path.exists() {
                eprintln!("Skipping test due to test data not found: {}\nDownload test data from: https://github.com/p-sira/ellip/tree/main/tests/data", $file_path);
                return;
            }

            let file = File::open($file_path).expect("Cannot open file");
            let reader = BufReader::new(file);

            let mut total_cases = 0;
            let mut test_fail = 0;
            let mut worst_line = 0;
            let mut worst_params: Vec<$t> = Vec::new();
            let mut worst_errors: [$t;5] = [0.0; 5];

            for (line_number, line) in reader.lines().enumerate() {
                let line = line.expect("Cannot read line");
                total_cases += 1;

                let parts: Vec<&str> = line.split_whitespace().collect();
                let inputs: Vec<$t> = parts[..parts.len() - 1]
                    .iter()
                    .map(|&v| $t::from_str(v).expect("Cannot parse input(s) as a number"))
                    .collect();

                let expected: $t = $t::from_str(parts[parts.len() - 1])
                    .expect("Cannot parse expected value as a number");

                let result = $func(&inputs);

                if result.is_nan() {
                    panic! ("Test failed on line {}: input = {:?}, expected = {:?}, got = NAN",
                    line_number + 1,
                    inputs,
                    expected,)
                }

                let error = (result - expected).abs();
                let tol = $t::from($atol) + $t::from($rtol) * expected.abs();
                if error > tol {
                    let rel = error / expected.abs();
                    test_fail += 1;
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

            if test_fail > 0 {
                panic!(
                    "Failed {}/{} cases. Worst on line {}: input = {:?}, expected = {:?}, got = {:?} \n error = {:?}, rel = {:?}, abs = {:?}, rtol = {:?}, atol = {:?}",
                    test_fail,
                    total_cases,
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
