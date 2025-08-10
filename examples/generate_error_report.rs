/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use std::fmt::Debug;

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

fn parse_wolfram_str<T: Float>(s: &str) -> Result<T, StrErr> {
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
    p99: f64,
}

impl Stats {
    pub fn nan() -> Self {
        return Stats {
            mean: f64::NAN,
            median: f64::NAN,
            variance: f64::NAN,
            max: f64::NAN,
            p99: f64::NAN,
        };
    }

    pub fn from_vec(v: &Vec<f64>) -> Self {
        let mut valids: Vec<f64> = v
            .iter()
            .filter(|x| !x.is_nan() && x.is_finite())
            .cloned()
            .collect();

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

        // Calculate P99 error
        let p99_pos = count as f64 * 0.99;
        let p99_pos_low = p99_pos as usize;
        let p99_frac = p99_pos - p99_pos_low as f64;
        let p99 = if p99_pos_low + 1 < count {
            valids[p99_pos_low] * (1.0 - p99_frac) + valids[p99_pos_low + 1] * p99_frac
        } else {
            valids[p99_pos_low]
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
            p99,
        }
    }
}

/// Reads wolfram data from file and returns a vector of Cases
fn read_wolfram_data<T: Float>(file_path: &str) -> Result<Vec<Case<T>>, StrErr> {
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

fn compute_errors_from_cases<T: Float + Debug>(
    func: &dyn Fn(&Vec<T>) -> T,
    cases: Vec<Case<T>>,
) -> Vec<f64> {
    cases
        .iter()
        .map(|case| {
            if case.expected.is_finite() {
                let res = func(&case.inputs);
                let err = rel_err(res, case.expected);
                // if err > 20.0 {
                //     println!(
                //         "Using parameters: {:?}, got={:?}, actual={:?} (error={:.2})",
                //         &case.inputs, res, case.expected, err
                //     );
                // }
                err
            } else {
                f64::NAN
            }
        })
        .collect()
}

fn format_float(value: &f64) -> String {
    if *value >= 1e5 {
        format!("{:.2e}", value)
    } else {
        format!("{:.2}", value)
    }
}

fn format_mu(value: &u64) -> String {
    if *value >= 10000 {
        format!("{:e}", value)
    } else {
        format!("{}", value)
    }
}

#[derive(Tabled)]
struct ErrorEntry<'a> {
    #[tabled(rename = "Function")]
    name: &'a str,
    #[tabled(rename = "Mean (ε)", display = "format_float")]
    mean: f64,
    #[tabled(rename = "Median (ε)", display = "format_float")]
    median: f64,
    #[tabled(rename = "P99 (ε)", display = "format_float")]
    p99: f64,
    #[tabled(rename = "Max (ε)", display = "format_float")]
    max: f64,
    #[tabled(rename = "Variance (ε²)", display = "format_float")]
    variance: f64,
    #[tabled(rename = "μ (ε²)", display = "format_mu")]
    mu: u64,
}

fn generate_error_entry_from_file<T: Float + Debug>(
    file_path: &str,
    func: &dyn Fn(&Vec<T>) -> T,
) -> Stats {
    let result = read_wolfram_data(file_path);
    match result {
        Ok(cases) => Stats::from_vec(&compute_errors_from_cases(func, cases)),
        Err(_) => Stats::nan(),
    }
}

fn generate_error_table(entries: &[(&str, Stats, u64)]) -> String {
    let rows: Vec<ErrorEntry> = entries
        .iter()
        .map(|(name, stats, mu)| ErrorEntry {
            name,
            mean: stats.mean,
            median: stats.median,
            p99: stats.p99,
            max: stats.max,
            variance: stats.variance,
            mu: *mu,
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
            let ans = $func(args[0], args[1], args[2]).unwrap_or_else(|e| {
                println!("{}: {}, {}, {}", e, args[0], args[1], args[2]);
                panic!()
            });
            if ans.is_nan() {
                println!("{}, {}, {}", args[0], args[1], args[2]);
                panic!("NAN")
            }
            ans
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
    ($file_name: expr, $name: expr, $func: expr, $arg_count: tt, $mu: expr) => {{
        func_wrapper!($func, $arg_count);

        let file_path = concat!["tests/data/", $file_name, ".csv"];

        (
            $name,
            generate_error_entry_from_file(&file_path, &wrapped_func),
            $mu,
        )
    }};
}

fn get_env() -> [String; 3] {
    use std::process::Command;

    let rust_version = {
        let output = Command::new("rustc")
            .arg("--version")
            .output()
            .unwrap()
            .stdout;
        String::from_utf8_lossy(&output)
            .split_whitespace()
            .collect::<Vec<&str>>()[1]
            .to_owned()
    };

    let platform = {
        let output = Command::new("rustc").arg("-vV").output().unwrap().stdout;
        String::from_utf8_lossy(&output)
            .lines()
            .find(|line| line.starts_with("host:"))
            .unwrap()
            .strip_prefix("host: ")
            .unwrap()
            .to_owned()
    };

    let ellip_version: String = {
        let output = Command::new("cargo")
            .args(["tree", "--invert", "--package", "ellip"])
            .output()
            .unwrap()
            .stdout;
        String::from_utf8_lossy(&output)
            .lines()
            .next()
            .and_then(|line| line.strip_prefix("ellip v"))
            .unwrap()
            .split_whitespace()
            .next()
            .unwrap()
            .to_owned()
    };

    [rust_version, platform, ellip_version]
}

fn main() {
    let [rust_version, platform, ellip_version] = get_env();

    let lines = [
         "# Testing",
         "This report presents the accuracy of the ellip crate using [**symmetric relative error**](https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/relative_error.html)",
         "metric. Errors are expressed in units of machine epsilon (ε).",
         "The test data spans the domain of each function up to **μ** to avoid approaching the function's limit.",
         "The reference values are computed using [**Wolfram Engine**](https://www.wolfram.com/engine/).",
         "You can find the scripts in the directory [tests/wolfram/](https://github.com/p-sira/ellip/blob/main/tests/wolfram/).",
         &format!("This report is generated on {} rustc {} using ellip v{} at `f64` precision (ε≈2.22e-16).", platform, rust_version, ellip_version),
         "",
         "## Legendre's Complete Elliptic Integrals",
         &generate_error_table(&[
             get_entry!("wolfram/ellipk_data", "ellipk", ellipk, 1, 1),
             get_entry!("wolfram/ellipk_neg", "ellipk (Neg m)", ellipk, 1, 1),
             get_entry!("wolfram/ellipe_data", "ellipe", ellipe, 1, 1),
             get_entry!("wolfram/ellipe_neg", "ellipe (Neg m)", ellipe, 1, 1),
             get_entry!("wolfram/ellippi_data", "ellippi", ellippi, 2, 1),
             get_entry!("wolfram/ellippi_neg", "ellippi (Neg m)", ellippi, 2, 1),
             get_entry!("wolfram/ellippi_pv", "ellippi (p.v.)", ellippi, 2, 50),
         ]),
         "",
         "## Legendre's Incomplete Elliptic Integrals",
         &generate_error_table(&[
             get_entry!("wolfram/ellipf_data", "ellipf", ellipf, 2, 1),
             get_entry!("wolfram/ellipf_neg", "ellipf (Neg m)", ellipf, 2, 1),
             get_entry!("wolfram/ellipeinc_data", "ellipeinc", ellipeinc, 2, 1),
             get_entry!("wolfram/ellipeinc_neg", "ellipeinc (Neg m)", ellipeinc, 2, 1),
             get_entry!("wolfram/ellippiinc_data", "ellippiinc", ellippiinc, 3, 1),
             get_entry!("wolfram/ellippiinc_neg", "ellippiinc (Neg m)", ellippiinc, 3, 1),
             get_entry!("wolfram/ellippiinc_pv", "ellippiinc (p.v.)", ellippiinc, 3, 1),
         ]),
         "",
         "## Bulirsch's Elliptic Integrals",
         "Bulirsh's elliptic integrals are not natively implemented in Wolfram Engine.",
         "Nevertheless, some of the integrals can be converted to their Legendre's counterpart,",
         "which are available on Wolfram Engine. However, for `cel` and `el2`, their values cannot",
         "be mapped simply. Hence, the reference values are generated using the functions",
         "submitted by Jan Mangaldan on [Wolfram Function Repository](https://resources.wolframcloud.com/FunctionRepository/).",
         "As for `cel2`, it is mapped to `cel` with p=1.",
         &generate_error_table(&[
             get_entry!("wolfram/cel_data", "cel", cel, 4, 1),
             get_entry!("wolfram/cel_pv", "cel (p.v.)", cel, 4, 1),
             get_entry!("wolfram/cel1_data", "cel1", cel1, 1, 1),
             get_entry!("wolfram/cel2_data", "cel2", cel2, 3, 1),
             get_entry!("wolfram/el1_data", "el1", el1, 2, 1),
             get_entry!("wolfram/el2_data", "el2", el2, 4, 1),
             get_entry!("wolfram/el3_data", "el3", el3, 3, 50),
             get_entry!("wolfram/el3_pv", "el3 (p.v.)", el3, 3, 50),
         ]),
         "",
         "## Carlson's Symmetric Elliptic Integrals",
         "",
         &generate_error_table(&[
             get_entry!("wolfram/elliprf_data", "elliprf", elliprf, 3, 1),
             get_entry!("wolfram/elliprg_data", "elliprg", elliprg, 3, 1),
             get_entry!("wolfram/elliprj_data", "elliprj", elliprj, 4, 50),
             get_entry!("wolfram/elliprj_pv", "elliprj (p.v.)", elliprj, 4, 10000000000),
             get_entry!("boost/elliprj_pv_small", "elliprj (p.v., small*)", elliprj, 4, 50),
             get_entry!("wolfram/elliprc_data", "elliprc", elliprc, 2, 1),
             get_entry!("wolfram/elliprc_pv", "elliprc (p.v.)", elliprc, 2, 1),
             get_entry!("wolfram/elliprd_data", "elliprd", elliprd, 3, 50),
         ]),
         "",
         "*small: Results compared with Boost Math implementation without promoting double, i.e, computed purely using `f64`.",
         "",
         "Current implementation of `elliprj` is less numerically stable in p.v. cases,",
         "as seen by large errors in the non-small test cases. That said, Ellip's results are consistent",
         "with Boost Math when limited to same precision (See [tests/data/boost/carlson.cpp](https://github.com/p-sira/ellip/blob/main/tests/data/boost/carlson.cpp)).",
         "Since the function is convergent, such errors can be mitigated when Rust's `f128` is released."
     ];

    use std::fs::File;
    use std::io::Write;
    let path = "tests/README.md";
    let mut file = File::create(path).unwrap();

    lines
        .iter()
        .for_each(|line| writeln!(file, "{}", line).expect("Cannot write line"));

    println!("Error report generated: {path}");
}
