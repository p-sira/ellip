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

    pub fn from_iter<I: Iterator<Item = f64>>(iter: I) -> Self {
        let mut count = 0;
        let mut sum = 0.0;
        let mut max = f64::NEG_INFINITY;
        let mut values = Vec::new();
        
        for x in iter {
            if !x.is_nan() && x.is_finite() {
                count += 1;
                sum += x;
                max = max.max(x);
                values.push(x);
            }
        }
        
        if count == 0 {
            return Self::nan();
        }
        
        let mean = sum / count as f64;
        
        values.sort_by(|a, b| a.partial_cmp(b).unwrap());
        
        // Calculate median
        let median = if count % 2 == 0 {
            (values[count / 2 - 1] + values[count / 2]) / 2.0
        } else {
            values[count / 2]
        };
        
        // Calculate variance
        let variance = values
            .iter()
            .map(|x| {
                let diff = mean - x;
                diff * diff
            })
            .sum::<f64>() / count as f64;
            
        Stats {
            mean,
            median,
            variance,
            max,
        }
    }

    #[allow(dead_code)]
    pub fn from_vec(v: &Vec<f64>) -> Self {
        Self::from_iter(v.iter().cloned())
    }
}

/// Reads wolfram data from file and returns a vector of Cases
fn read_wolfram_data<T: Float + 'static>(
    file_path: &str,
) -> Box<dyn Iterator<Item = Case<T>>> {
    use csv::ReaderBuilder;
    use std::fs::File;
    use std::path::Path;
    use std::io::BufReader;

    let path = Path::new(file_path);
    let file = match File::open(path) {
        Ok(f) => f,
        Err(_) => {
            eprintln!(
                "Test data not found: {file_path}\nDownload from: https://github.com/p-sira/ellip/tree/main/tests/data"
            );
            return Box::new(std::iter::empty());
        }
    };

    let reader = ReaderBuilder::new()
        .has_headers(false)
        .from_reader(BufReader::new(file));

    let iter = reader
        .into_records()
        .filter_map(|result| match result {
            Ok(row) => {
                let len = row.len();
                let inputs = row
                    .iter()
                    .take(len - 1)
                    .map(parse_wolfram_str)
                    .collect::<Result<Vec<T>, _>>()
                    .ok()?;
                let expected = parse_wolfram_str(&row[len - 1]).ok()?;
                Some(Case { inputs, expected })
            }
            Err(_) => None,
        });

    Box::new(iter)
}

fn compute_errors_from_cases<T: Float + 'static, I: Iterator<Item = Case<T>>>(
    func: &'static dyn Fn(&Vec<T>) -> T,
    cases: I,
) -> impl Iterator<Item = f64> {
    cases.map(move |case| {
        if case.expected.is_finite() {
            let res = func(&case.inputs);
            rel_err(res, case.expected)
        } else {
            f64::NAN
        }
    })
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

fn generate_error_entry_from_file<T: Float + 'static>(file_path: &str, func: &'static dyn Fn(&Vec<T>) -> T) -> Stats {
    let cases = read_wolfram_data(file_path);
    let errors = compute_errors_from_cases(func, cases);
    Stats::from_iter(errors)
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
        "Since the test relied on fixed `f64` precision, the precision deteriorates at extreme values (too big or too small).",
        "However, the calculation via Wolfram Engine utilized arbitary precision up to 50 digits.",
        "This could lead to decrepancies in the results.",
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
            get_entry!("ellippiinc_data", "ellippiinc", ellippiinc, 3),
            get_entry!("ellippiinc_neg", "ellippiinc (Negative m)", ellippiinc, 3),
            get_entry!("ellippiinc_limit", "ellippiinc (Near limits)", ellippiinc, 3),
        ]),
        "",
        "## Bulirsch's Elliptic Integrals",
        "For Bulirsh's elliptic integrals, the reference values are generated using the function",
        "submitted by Jan Mangaldan on [Wolfram Function Repository](https://resources.wolframcloud.com/FunctionRepository/).",

        &generate_error_table(&[
            get_entry!("cel_data", "cel", cel, 4),
            get_entry!("cel_pv", "cel (p.v.)", cel, 4),
            get_entry!("el1_data", "el1", el1, 2),
            get_entry!("el1_limit", "el1 (Near limits)", el1, 2),
            get_entry!("el2_data", "el2", el2, 4),
            get_entry!("el2_limit", "el2 (Near limits)", el2, 4),
            get_entry!("el3_data", "el3", el3, 3),
            get_entry!("el3_pv", "el3 (p.v.)", el3, 3),
            get_entry!("el3_limit", "el3 (Near limits)", el3, 3),
        ]),
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
