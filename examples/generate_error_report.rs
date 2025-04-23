/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::*;
use num_traits::Float;
use tabled::{settings::Style, Table, Tabled};

/// Calculates max relative error in unit of epsilon
fn rel_err<T: Float>(a: T, b: T) -> f64 {
    let err_a = ((a - b) / a).abs();
    let err_b = ((a - b) / b).abs();

    let rel_err = err_a.max(err_b);
    (rel_err / T::epsilon())
        .to_f64()
        .expect("Cannot convert to f64")
}

fn parse_wolfram_str<T: Float>(s: &str) -> Result<T, &'static str> {
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

#[derive(Debug, Clone)]
struct Case<T: Float> {
    inputs: Vec<T>,
    expected: T,
}

struct Stats {
    mean: f64,
    median: f64,
    variance: f64,
    max: f64,
}

impl Stats {
    pub fn nan() -> Self {
        return Stats {
            mean: f64::NAN,
            median: f64::NAN,
            variance: f64::NAN,
            max: f64::NAN,
        };
    }

    pub fn from_vec(v: &Vec<f64>) -> Self {
        let mut valids: Vec<f64> = v.iter().filter(|x| x.is_normal()).cloned().collect();

        if valids.is_empty() {
            return Self::nan();
        }

        valids.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let sum: f64 = valids.iter().sum();
        let count = valids.len();
        let mean = sum / count as f64;

        // Calculate median
        let median = if count % 2 == 0 {
            (valids[count / 2 - 1] + valids[count / 2]) / 2.0
        } else {
            valids[count / 2]
        };

        // Calculate variance
        let variance = valids
            .iter()
            .map(|x| {
                let diff = mean - x;
                diff * diff
            })
            .sum::<f64>()
            / count as f64;

        let max = *valids
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(&f64::NAN);

        Stats {
            mean,
            median,
            variance,
            max,
        }
    }
}

/// Reads wolfram data from file and returns a vector of Cases
fn read_wolfram_data<T: Float>(file_path: &str) -> Result<Vec<Case<T>>, &'static str> {
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

fn compute_errors_from_cases<T: Float>(
    func: &dyn Fn(&Vec<T>) -> T,
    cases: Vec<Case<T>>,
) -> Vec<f64> {
    cases
        .iter()
        .map(|case| {
            if case.expected.is_finite() {
                let res = func(&case.inputs);
                rel_err(res, case.expected)
            } else {
                f64::NAN
            }
        })
        .collect()
}

fn format_float(value: &f64) -> String {
    format!("{:.2}", value)
}

#[derive(Tabled)]
struct ErrorEntry<'a> {
    #[tabled(rename = "Function")]
    name: &'a str,
    #[tabled(rename = "Mean (ε)", display = "format_float")]
    mean: f64,
    #[tabled(rename = "Median (ε)", display = "format_float")]
    median: f64,
    #[tabled(rename = "Variance (ε²)", display = "format_float")]
    variance: f64,
    #[tabled(rename = "Max (ε)", display = "format_float")]
    max: f64,
}

fn generate_error_entry_from_file<T: Float>(file_path: &str, func: &dyn Fn(&Vec<T>) -> T) -> Stats {
    let result = read_wolfram_data(file_path);
    match result {
        Ok(cases) => Stats::from_vec(&compute_errors_from_cases(func, cases)),
        Err(_) => Stats::nan(),
    }
}

fn generate_error_table(entries: &[(&str, Stats)]) -> String {
    let rows: Vec<ErrorEntry> = entries
        .iter()
        .map(|(name, stats)| ErrorEntry {
            name,
            mean: stats.mean,
            median: stats.median,
            variance: stats.variance,
            max: stats.max,
        })
        .collect();

    Table::new(rows).with(Style::markdown()).to_string()
}

macro_rules! func_wrapper {
    ($func:expr, 1) => {
        fn wrapped_func(args: &Vec<f64>) -> f64 {
            $func(args[0]).unwrap()
        }
    };
    ($func:expr, 2) => {
        fn wrapped_func(args: &Vec<f64>) -> f64 {
            $func(args[0], args[1]).unwrap()
        }
    };
    ($func:expr, 3) => {
        fn wrapped_func(args: &Vec<f64>) -> f64 {
            $func(args[0], args[1], args[2]).unwrap()
        }
    };
    ($func:expr, 4) => {
        fn wrapped_func(args: &Vec<f64>) -> f64 {
            $func(args[0], args[1], args[2], args[3]).unwrap()
        }
    };
    ($func:expr, $_:expr) => {
        fn wrapped_func(_args: &Vec<f64>) -> f64 {
            panic!("Unsupported number of arguments")
        }
    };
}

macro_rules! get_entry {
    ($file_name: expr, $name: expr, $func: expr, $arg_count: tt) => {{
        func_wrapper!($func, $arg_count);

        let file_path = ["./tests/data/wolfram/", $file_name, ".csv"].concat();

        (
            $name,
            generate_error_entry_from_file(&file_path, &wrapped_func),
        )
    }};
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    use std::process::Command;

    use std::fs::File;
    use std::io::Write;

    let rust_version = {
        let output = Command::new("rustc").arg("--version").output()?.stdout;
        String::from_utf8_lossy(&output)
            .split_whitespace()
            .collect::<Vec<&str>>()[1]
            .to_owned()
    };

    let platform = {
        let output = Command::new("rustc").arg("-vV").output()?.stdout;
        String::from_utf8_lossy(&output)
            .lines()
            .find(|line| line.starts_with("host:"))
            .unwrap()
            .strip_prefix("host: ")
            .unwrap()
            .to_owned()
    };

    let ellip_version = env!("CARGO_PKG_VERSION");

    let lines = [
        "# Testing",
        "This report presents the accuracy of the ellip crate using [symmetric relative error](https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/relative_error.html)",
        "metric. The errors are expressed in units of machine epsilon (ε).",
        "The reference values are computed using [Wolfram Engine](https://www.wolfram.com/engine/),",
        "You can find the scripts in the directory [tests/wolfram/](https://github.com/p-sira/ellip/blob/main/tests/wolfram/).",
        &format!("This report is generated on {} rustc {} using ellip v{} at `f64` precision (ε≈2.22e-16).", platform, rust_version, ellip_version),
        "",
        "## Legendre's Complete Elliptic Integrals",
        &generate_error_table(&[
            get_entry!("ellipk_data", "ellipk", ellipk, 1),
            get_entry!("ellipk_neg", "ellipk (Negative m)", ellipk, 1),
            get_entry!("ellipk_limit", "ellipk (Near limits)", ellipk, 1),
            get_entry!("ellipe_data", "ellipe", ellipe, 1),
            get_entry!("ellipe_neg", "ellipe (Negative m)", ellipe, 1),
            get_entry!("ellipe_limit", "ellipe (Near limits)", ellipe, 1),
            get_entry!("ellippi_data", "ellippi", ellippi, 2),
            get_entry!("ellippi_neg", "ellippi (Negative m)", ellippi, 2),
            get_entry!("ellippi_pv", "ellippi (p.v.)", ellippi, 2),
            get_entry!("ellippi_limit", "ellippi (Near limits)", ellippi, 2),
        ]),
        "",
        "## Legendre's Incomplete Elliptic Integrals",
        &generate_error_table(&[
            get_entry!("ellipf_data", "ellipf", ellipf, 2),
            get_entry!("ellipf_neg", "ellipf (Negative m)", ellipf, 2),
            get_entry!("ellipf_limit", "ellipf (Near limits)", ellipf, 2),
            get_entry!("ellipeinc_data", "ellipeinc", ellipeinc, 2),
            get_entry!("ellipeinc_neg", "ellipeinc (Negative m)", ellipeinc, 2),
            get_entry!("ellipeinc_limit", "ellipeinc (Near limits)", ellipeinc, 2),
        ]),
        "",
        "## Bulirsch's Complete Elliptic Integrals",
        "",
        "## Bulirsch's Incomplete Elliptic Integrals",
        "",
        "## Carlson's Symmetric Elliptic Integrals",
        "",
    ];

    let path = "tests/README.md";
    let mut file = File::create(path)?;

    lines
        .iter()
        .for_each(|line| writeln!(file, "{}", line).expect("Cannot write line"));

    println!("Error report generated: {path}");

    Ok(())
}
