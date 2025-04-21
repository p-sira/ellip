/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use std::fmt::{Display, Error, Formatter};

use num_traits::Float;
use serde::Deserialize;

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

/*
 * WolframFloat
*/
#[derive(Debug)]
pub enum WolframFloat {
    Value(f64),
    Infinity,
    NegativeInfinity,
    ComplexInfinity,
    Indeterminate,
    NaN,
}

#[allow(dead_code)]
impl WolframFloat {
    pub fn to_float<T: Float>(&self) -> T {
        match *self {
            WolframFloat::Value(v) => T::from(v).unwrap_or_else(|| nan!()),
            WolframFloat::Infinity => inf!(),
            WolframFloat::NegativeInfinity => neg_inf!(),
            WolframFloat::ComplexInfinity => inf!(),
            WolframFloat::Indeterminate => nan!(),
            WolframFloat::NaN => nan!(),
        }
    }
}

impl<'de> Deserialize<'de> for WolframFloat {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s: &str = Deserialize::deserialize(deserializer)?;
        match s {
            "Infinity" => Ok(WolframFloat::Infinity),
            "-Infinity" => Ok(WolframFloat::NegativeInfinity),
            "ComplexInfinity" => Ok(WolframFloat::ComplexInfinity),
            "Indeterminate" => Ok(WolframFloat::Indeterminate),
            "NaN" => Ok(WolframFloat::NaN),
            _ => s
                .parse::<f64>()
                .map(WolframFloat::Value)
                .map_err(serde::de::Error::custom),
        }
    }
}

impl Display for WolframFloat {
    fn fmt(&self, f: &mut Formatter) -> Result<(), Error> {
        match *self {
            WolframFloat::Value(v) => write!(f, "{}", v),
            WolframFloat::Infinity => write!(f, "Infinity"),
            WolframFloat::NegativeInfinity => write!(f, "NegativeInfinity"),
            WolframFloat::ComplexInfinity => write!(f, "ComplexInfinity"),
            WolframFloat::Indeterminate => write!(f, "Indeterminate"),
            WolframFloat::NaN => write!(f, "NaN"),
        }
    }
}
