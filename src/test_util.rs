/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use std::{
    fmt::{Display, Error, Formatter},
    str::FromStr,
};

use num_traits::Float;

#[macro_export]
macro_rules! compare_test_data {
    ($func: expr, $inputs :expr, $expected: expr, $t: ident, $rtol:expr, $atol:expr) => {
        let mut test_fail = 0;
        let mut worst_line = 0;
        let mut worst_params: Vec<$t> = Vec::new();
        let mut worst_errors: [$t;5] = [0.0; 5];

        for (i, (inp, exp)) in $inputs.iter().zip($expected).enumerate() {
            let result = $func(&inp);

            if exp.is_nan() && result.is_nan() {
                continue;
            }

            if result.is_nan() {
                eprintln! ("Test failed on case {}: input = {:?}, expected = {:?}, got = NAN",
                    i + 1,
                    inp,
                    exp,);
                test_fail += 1;
            }

            if !result.is_finite() {
                if result != exp {
                    eprintln! ("Test failed on case {}: input = {:?}, expected = {:?}, got = {:?}",
                        i + 1,
                        inp,
                        exp,
                        result);
                    test_fail += 1;
                }
            }

            let error = (result - exp).abs();
            let tol = $t::from($atol) + $t::from($rtol) * exp.abs();
            if error > tol {
                let rel = error / exp.abs();
                test_fail += 1;
                if rel > worst_errors[1] {
                    worst_line = i + 1;
                    worst_errors = [error, rel, tol, exp, result];
                    worst_params = inp.clone();
                }
                eprintln!(
                    "Test failed on case {}: input = {:?}, expected = {:?}, got = {:?} \n error = {:?}, rel = {:?}, abs = {:?}, rtol = {:?}, atol = {:?}\n",
                    i + 1,
                    inp,
                    exp,
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
                "Failed {}/{} cases. Worst on case {}: input = {:?}, expected = {:?}, got = {:?} \n error = {:?}, rel = {:?}, abs = {:?}, rtol = {:?}, atol = {:?}",
                test_fail,
                $inputs.len(),
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
    };
}

#[macro_export]
macro_rules! compare_test_data_boost {
    ($filename:expr, $func:expr, $rtol:expr) => {
        compare_test_data_boost!($filename, $func, f64, $rtol, 0.0)
    };
    ($filename:expr, $func:expr, $t: ident, $rtol:expr) => {
        compare_test_data_boost!($filename, $func, f64, $rtol, 0.0)
    };
    ($filename:expr, $func:expr, $rtol:expr, $atol:expr) => {
        compare_test_data_boost!($filename, $func, f64, $rtol, $atol)
    };
    ($filename:expr, $func:expr, $t: ident, $rtol:expr, $atol:expr) => {
        {
            use crate::compare_test_data;

            use std::fs::File;
            use std::io::{BufRead, BufReader};
            use std::path::Path;
            use std::str::FromStr;

            let path = Path::new("./tests/data/boost").join($filename);
            if !path.exists() {
                eprintln!(
                    "Skipping test due to test data not found: {}\nDownload test data from: https://github.com/p-sira/ellip/tree/main/tests/data",
                    path.to_str().unwrap_or("{path not utf-8}")
                );
                return;
            }

            let file = File::open(path).expect("Cannot open file");
            let reader = BufReader::new(file);

            let mut inputs: Vec<Vec<$t>> = Vec::new();
            let mut expected: Vec<$t> = Vec::new();

            reader.lines().for_each(|line| {
                let line = line.expect("Cannot read line");
                let parts: Vec<&str> = line.split_whitespace().collect();
                inputs.push(
                    parts[..parts.len() - 1]
                        .iter()
                        .map(|&v| $t::from_str(v).expect("Cannot parse input(s) as a number"))
                        .collect(),
                );
                expected.push(
                    $t::from_str(parts[parts.len() - 1]).expect("Cannot parse expected value as a number"),
                );
            });

            compare_test_data!($func, inputs, expected, $t, $rtol, $atol);
        }
    };
}

#[macro_export]
macro_rules! compare_test_data_wolfram {
    ($filename:expr, $func:expr, $rtol:expr) => {
        compare_test_data_wolfram!($filename, $func, f64, $rtol, 0.0)
    };
    ($filename:expr, $func:expr, $t: ident, $rtol:expr, $atol:expr) => {{
        {
            use crate::compare_test_data;
            use crate::test_util::WolframFloat;
            use std::str::FromStr;

            use std::fs::File;
            use std::path::Path;
            use csv::ReaderBuilder;

            let path = Path::new("./tests/data/wolfram").join($filename);
            if !path.exists() {
                eprintln!(
                    "Skipping test due to test data not found: {}\nDownload test data from: https://github.com/p-sira/ellip/tree/main/tests/data",
                    path.to_str().unwrap_or("{path not utf-8}")
                );
                return;
            }

            let file = File::open(path).expect("Cannot open file");
            let mut reader = ReaderBuilder::new().has_headers(false).from_reader(file);

            let mut inputs: Vec<Vec<$t>> = Vec::new();
            let mut expected: Vec<$t> = Vec::new();

            for record in reader.records() {
                let row = record.expect("Cannot read row");
                let len = row.len();

                inputs.push(
                    row.iter().take(len - 1).map(|s| {
                        WolframFloat::from_str(s).expect("Cannot parse input(s) as a number")
                        .to_float::<$t>()
                    }).collect()
                );
                expected.push(
                    WolframFloat::from_str(&row[len - 1]).expect("Cannot parse expected value as a number")
                    .to_float::<$t>(),
                );
            }

            compare_test_data!($func, inputs, expected, $t, $rtol, $atol);
        }
    }};
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
#[allow(dead_code)]
pub enum WolframFloat {
    Float64(f64),
    Float32(f32),
    Infinity,
    NegativeInfinity,
    ComplexInfinity,
    Indeterminate,
    NaN,
}

impl WolframFloat {
    pub fn to_float<T: Float>(&self) -> T {
        match *self {
            WolframFloat::Float64(v) => T::from(v).unwrap_or_else(|| nan!()),
            WolframFloat::Float32(v) => T::from(v).unwrap_or_else(|| nan!()),
            WolframFloat::Infinity => inf!(),
            WolframFloat::NegativeInfinity => neg_inf!(),
            WolframFloat::ComplexInfinity => inf!(),
            WolframFloat::Indeterminate => nan!(),
            WolframFloat::NaN => nan!(),
        }
    }
}

impl FromStr for WolframFloat {
    type Err = std::num::ParseFloatError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim() {
            "Infinity" => Ok(WolframFloat::Infinity),
            "-Infinity" => Ok(WolframFloat::NegativeInfinity),
            "ComplexInfinity" => Ok(WolframFloat::ComplexInfinity),
            "Indeterminate" => Ok(WolframFloat::Indeterminate),
            "NaN" => Ok(WolframFloat::NaN),
            _ => Ok(WolframFloat::Float64(s.parse::<f64>()?)),
        }
    }
}

impl Display for WolframFloat {
    fn fmt(&self, f: &mut Formatter) -> Result<(), Error> {
        match *self {
            WolframFloat::Float64(v) => write!(f, "{}", v),
            WolframFloat::Float32(v) => write!(f, "{}", v),
            WolframFloat::Infinity => write!(f, "Infinity"),
            WolframFloat::NegativeInfinity => write!(f, "NegativeInfinity"),
            WolframFloat::ComplexInfinity => write!(f, "ComplexInfinity"),
            WolframFloat::Indeterminate => write!(f, "Indeterminate"),
            WolframFloat::NaN => write!(f, "NaN"),
        }
    }
}
