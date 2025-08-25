/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use num_traits::Float;
use std::fmt::Debug;
use tabled::{Table, Tabled, settings::Style};

use crate::stats::Stats;

/// Calculates max relative error in unit of epsilon
pub fn rel_err<T: Float>(a: T, b: T) -> f64 {
    let err_a = ((a - b) / a).abs();
    let err_b = ((a - b) / b).abs();

    let rel_err = err_a.max(err_b);
    (rel_err / T::epsilon())
        .to_f64()
        .expect("Cannot convert to f64")
}

#[derive(Debug, Clone)]
pub struct Case<T: Float> {
    pub inputs: Vec<T>,
    pub expected: T,
}

pub fn compute_errors_from_cases<T: Float + Debug>(
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
    if value.is_nan() {
        "NAN".to_string()
    } else if *value >= 1e3 {
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
pub struct ErrorEntry<'a> {
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

pub fn generate_error_entry_from_file<T: Float + Debug>(
    file_path: &str,
    func: &dyn Fn(&Vec<T>) -> T,
) -> Stats {
    let result = crate::parser::read_wolfram_data(file_path);
    match result {
        Ok(cases) => Stats::from_vec(&compute_errors_from_cases(func, cases)),
        Err(_) => Stats::nan(),
    }
}

pub fn generate_error_table(entries: &[(&str, Stats, u64)]) -> String {
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

#[macro_export]
macro_rules! get_entry {
    ($file_name: expr, $name: expr, $func: expr, $arg_count: tt, $mu: expr) => {{
        ellip_dev_utils::func_wrapper!($func, $arg_count);

        let file_path = concat!["tests/data/", $file_name, ".csv"];

        (
            $name,
            ellip_dev_utils::test_report::generate_error_entry_from_file(&file_path, &wrapped_func),
            $mu,
        )
    }};
}

fn format_performance(value: &f64) -> String {
    if value.is_nan() {
        "NAN".to_string()
    } else if *value < 1000.0 {
        format!("{:.1} ns", value)
    } else if *value < 1_000_000.0 {
        format!("{:.1} μs", value / 1000.0)
    } else {
        format!("{:.1} ms", value / 1_000_000.0)
    }
}

#[derive(Tabled)]
pub struct SummaryEntry<'a> {
    #[tabled(rename = "Function")]
    name: &'a str,
    #[tabled(rename = "Median Error (ε)", display = "format_float")]
    median_error: f64,
    #[tabled(rename = "Max Error (ε)", display = "format_float")]
    max_error: f64,
    #[tabled(rename = "Mean Performance", display = "format_performance")]
    mean_performance: f64,
}

pub fn generate_summary_table(entries: &[(&str, Stats, f64)]) -> String {
    let rows: Vec<SummaryEntry> = entries
        .iter()
        .map(|(name, stats, perf)| SummaryEntry {
            name,
            median_error: stats.median,
            max_error: stats.max,
            mean_performance: *perf / stats.n as f64,
        })
        .collect();

    Table::new(rows).with(Style::markdown()).to_string()
}

#[macro_export]
macro_rules! get_summary_entry {
    ($group:expr, $name:expr, $func:expr, $arg_count:tt, $test_file_name:expr) => {{
        use ellip_dev_utils::{
            benchmark, file, parser, stats,
            test_report::{self, Case},
        };
        use std::path::Path;

        ellip_dev_utils::func_wrapper!($func, $arg_count);

        let test_paths = file::find_test_files($test_file_name, "wolfram");
        let cases = test_paths
            .iter()
            .flat_map(|test_path| parser::read_wolfram_data(test_path.to_str().unwrap()).unwrap())
            .collect::<Vec<Case<f64>>>();
        let stats = stats::Stats::from_vec(&test_report::compute_errors_from_cases(
            &wrapped_func,
            cases,
        ));

        // Criterion directory structure: target/criterion/<group>/<func>/new/estimates.json
        let criterion_path_buf = Path::new("target/criterion")
            .join($group)
            .join(stringify!($func))
            .join("new");
        let perf = benchmark::extract_criterion_mean(criterion_path_buf.to_str().unwrap_or(""))
            .unwrap_or(f64::NAN);

        ($name, stats, perf)
    }};
    ($group:expr, $name:expr, $func:expr, $arg_count:tt) => {{
        get_summary_entry! {$group, $name, $func, $arg_count, stringify!($func)}
    }};
}
