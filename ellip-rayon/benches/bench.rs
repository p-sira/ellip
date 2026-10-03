/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use std::io::Write;
use std::path::{Path, PathBuf};

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use ellip_dev_utils::{benchmark::extract_criterion_mean, parser, test_report::Case};
use itertools::izip;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};
use tabled::{Table, Tabled, settings::Style};

const MAX_THRESHOLD: usize = 10000;
const ACCEPT_THRESHOLD: f64 = 0.1;
const STEP: usize = 100;
const FILE_SAMPLE_SIZE: usize = 500;
const MAX_ITER: usize = 15;

pub fn linspace(start: usize, end: usize, n: usize) -> Vec<usize> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![start];
    }

    let delta = end - start - 1;
    (0..=(n - 1)).map(|i| start + i * delta / (n - 1)).collect()
}

fn get_cases(test_paths: &[PathBuf], n: usize) -> Vec<Case<f64>> {
    let cases = test_paths
        .iter()
        .flat_map(|test_path| parser::read_wolfram_data(test_path.to_str().unwrap()).unwrap())
        .collect::<Vec<Case<f64>>>();
    let base_cases: Vec<Case<f64>> = linspace(0, cases.len(), FILE_SAMPLE_SIZE)
        .iter()
        .map(|&i| cases[i].clone())
        .collect();
    base_cases.into_iter().cycle().take(n).collect()
}

fn extract_params(cases: &[Case<f64>], param_index: usize) -> Vec<f64> {
    cases.iter().map(|case| case.inputs[param_index]).collect()
}

pub fn find_test_files(function_name: &str, dir: &str) -> Vec<PathBuf> {
    let pattern = format!("../tests/data/{dir}/{function_name}_*.csv");
    let mut files: Vec<PathBuf> = glob::glob(&pattern)
        .unwrap()
        .map(|path| path.unwrap())
        .collect();
    if files.is_empty() {
        let pattern_root = format!("tests/data/{dir}/{function_name}_*.csv");
        files = glob::glob(&pattern_root)
            .unwrap()
            .map(|path| path.unwrap())
            .collect();
    }
    files
}

macro_rules! find_test_files {
    ($func:ident) => {
        find_test_files(stringify!($func), "wolfram")
    };
    ($func:ident : $test_file_name:expr) => {
        find_test_files($test_file_name, "wolfram")
    };
}

enum ExitCond {
    None,
    Converged,
    MinStep,
    MaxIter,
    MaxThres,
}

impl ExitCond {
    fn to_string(&self) -> String {
        match self {
            Self::None => "-".to_string(),
            Self::Converged => "Converged".to_string(),
            Self::MinStep => "Min step size".to_string(),
            Self::MaxIter => "Max iteration".to_string(),
            Self::MaxThres => "Max threshold".to_string(),
        }
    }
}

fn format_float(val: &f64) -> String {
    if val.is_nan() {
        "-".to_string()
    } else {
        format!("{:.3}", val)
    }
}

fn format_performance(val: &f64) -> String {
    let val = *val / 1e9;
    if val.is_nan() {
        "-".to_string()
    } else if val < 1e-6 {
        format!("{:.3} ns", val * 1e9)
    } else if val < 1e-3 {
        format!("{:.3} µs", val * 1e6)
    } else if val < 1.0 {
        format!("{:.3} ms", val * 1e3)
    } else {
        format!("{:.3} s", val)
    }
}

#[derive(Tabled)]
struct Record {
    #[tabled(rename = "Function")]
    function_name: String,
    #[tabled(rename = "Threshold")]
    threshold: usize,
    #[tabled(rename = "Parallel Time", display = "format_performance")]
    par_time: f64,
    #[tabled(rename = "Serial Time", display = "format_performance")]
    ser_time: f64,
    #[tabled(rename = "Ratio", display = "format_float")]
    ratio: f64,
    #[tabled(rename = "Exit Condition", display = "ExitCond::to_string")]
    exit_cond: ExitCond,
}

fn parse_performance(s: &str) -> f64 {
    let s = s.trim();
    if s == "-" {
        return f64::NAN;
    }
    if let Some(num) = s.strip_suffix("ns") {
        num.trim().parse::<f64>().unwrap_or(f64::NAN)
    } else if let Some(num) = s.strip_suffix("µs").or_else(|| s.strip_suffix("us")) {
        num.trim()
            .parse::<f64>()
            .map(|v| v * 1e3)
            .unwrap_or(f64::NAN)
    } else if let Some(num) = s.strip_suffix("ms") {
        num.trim()
            .parse::<f64>()
            .map(|v| v * 1e6)
            .unwrap_or(f64::NAN)
    } else if let Some(num) = s.strip_suffix("s") {
        num.trim()
            .parse::<f64>()
            .map(|v| v * 1e9)
            .unwrap_or(f64::NAN)
    } else {
        f64::NAN
    }
}

fn load_results(path: &Path) -> Vec<Record> {
    let Ok(content) = std::fs::read_to_string(path) else {
        return Vec::new();
    };
    let mut records = Vec::new();
    for line in content.lines().skip(2) {
        let line = line.trim();
        if line.is_empty() || !line.starts_with('|') {
            continue;
        }
        let parts: Vec<&str> = line.split('|').map(|s| s.trim()).collect();
        if parts.len() >= 7 {
            let function_name = parts[1].to_string();
            let Ok(threshold) = parts[2].parse::<usize>() else {
                continue;
            };
            let par_time = parse_performance(parts[3]);
            let ser_time = parse_performance(parts[4]);
            let ratio = parts[5].parse::<f64>().unwrap_or(f64::NAN);
            let exit_cond = match parts[6] {
                "Converged" => ExitCond::Converged,
                "Min step size" => ExitCond::MinStep,
                "Max iteration" => ExitCond::MaxIter,
                "Max threshold" => ExitCond::MaxThres,
                _ => ExitCond::None,
            };
            records.push(Record {
                function_name,
                threshold,
                par_time,
                ser_time,
                ratio,
                exit_cond,
            });
        }
    }
    records
}

fn update_or_push_record(results: &mut Vec<Record>, record: Record) {
    if let Some(pos) = results
        .iter()
        .position(|r| r.function_name == record.function_name)
    {
        results[pos] = record;
    } else {
        results.push(record);
    }
}

fn save_results(path: &Path, results: &[Record]) {
    let result_str = Table::new(results).with(Style::markdown()).to_string();
    let mut file = std::fs::File::create(path).unwrap();
    file.write_all(result_str.as_bytes()).unwrap();
}

fn get_criterion_root() -> PathBuf {
    if Path::new("target/criterion").exists() || !Path::new("../target/criterion").exists() {
        PathBuf::from("target/criterion")
    } else {
        PathBuf::from("../target/criterion")
    }
}

fn get_record_file() -> PathBuf {
    if Path::new("ellip-rayon/benches").exists() {
        PathBuf::from("ellip-rayon/benches/par_threshold.md")
    } else {
        PathBuf::from("benches/par_threshold.md")
    }
}

macro_rules! par_zip {
    ($a:expr) => {
        $a.par_iter()
    };
    ($a:expr, $b:expr) => {
        $a.par_iter().zip($b.par_iter())
    };
    ($a:expr, $b:expr, $c:expr) => {
        $a.par_iter()
            .zip($b.par_iter())
            .zip($c.par_iter())
            .map(|((a, b), c)| (a, b, c))
    };
    ($a:expr, $b:expr, $c:expr, $d:expr) => {
        $a.par_iter()
            .zip($b.par_iter())
            .zip($c.par_iter())
            .zip($d.par_iter())
            .map(|(((a, b), c), d)| (a, b, c, d))
    };
}

macro_rules! bench_func {
    (@run $group:expr, $func:ident, $test_paths:expr, $size:expr, 1) => {{
        let cases = get_cases($test_paths, $size);
        let a = extract_params(&cases, 0);
        $group.bench_with_input(BenchmarkId::new("serial", $size), &$size, |b, _| {
            b.iter(|| {
                let ans: Vec<_> = a.iter().map(|&x| ellip::$func(x).unwrap()).collect();
                std::hint::black_box(&ans);
            });
        });
        $group.bench_with_input(BenchmarkId::new("batch", $size), &$size, |b, _| {
            b.iter(|| {
                let ans: Vec<_> = a.par_iter().map(|&x| ellip::$func(x).unwrap()).collect();
                std::hint::black_box(&ans);
            });
        });
    }};
    (@run $group:expr, $func:ident, $test_paths:expr, $size:expr, 2) => {{
        let cases = get_cases($test_paths, $size);
        let a = extract_params(&cases, 0);
        let b = extract_params(&cases, 1);
        $group.bench_with_input(BenchmarkId::new("serial", $size), &$size, |bench, _| {
            bench.iter(|| {
                let ans: Vec<_> = izip!(&a, &b)
                    .map(|(&x, &y)| ellip::$func(x, y).unwrap())
                    .collect();
                std::hint::black_box(&ans);
            });
        });
        $group.bench_with_input(BenchmarkId::new("batch", $size), &$size, |bench, _| {
            bench.iter(|| {
                let ans: Vec<_> = par_zip!(&a, &b)
                    .map(|(&x, &y)| ellip::$func(x, y).unwrap())
                    .collect();
                std::hint::black_box(&ans);
            });
        });
    }};
    (@run $group:expr, $func:ident, $test_paths:expr, $size:expr, 3) => {{
        let cases = get_cases($test_paths, $size);
        let a = extract_params(&cases, 0);
        let b = extract_params(&cases, 1);
        let c = extract_params(&cases, 2);
        $group.bench_with_input(BenchmarkId::new("serial", $size), &$size, |bench, _| {
            bench.iter(|| {
                let ans: Vec<_> = izip!(&a, &b, &c)
                    .map(|(&x, &y, &z)| ellip::$func(x, y, z).unwrap())
                    .collect();
                std::hint::black_box(&ans);
            });
        });
        $group.bench_with_input(BenchmarkId::new("batch", $size), &$size, |bench, _| {
            bench.iter(|| {
                let ans: Vec<_> = par_zip!(&a, &b, &c)
                    .map(|(&x, &y, &z)| ellip::$func(x, y, z).unwrap())
                    .collect();
                std::hint::black_box(&ans);
            });
        });
    }};
    (@run $group:expr, $func:ident, $test_paths:expr, $size:expr, 4) => {{
        let cases = get_cases($test_paths, $size);
        let a = extract_params(&cases, 0);
        let b = extract_params(&cases, 1);
        let c = extract_params(&cases, 2);
        let d = extract_params(&cases, 3);
        $group.bench_with_input(BenchmarkId::new("serial", $size), &$size, |bench, _| {
            bench.iter(|| {
                let ans: Vec<_> = izip!(&a, &b, &c, &d)
                    .map(|(&x, &y, &z, &p)| ellip::$func(x, y, z, p).unwrap())
                    .collect();
                std::hint::black_box(&ans);
            });
        });
        $group.bench_with_input(BenchmarkId::new("batch", $size), &$size, |bench, _| {
            bench.iter(|| {
                let ans: Vec<_> = par_zip!(&a, &b, &c, &d)
                    .map(|(&x, &y, &z, &p)| ellip::$func(x, y, z, p).unwrap())
                    .collect();
                std::hint::black_box(&ans);
            });
        });
    }};
}

macro_rules! find_threshold {
    ($c:ident, $results:ident, $func:ident, $n_args:tt, $init_step:expr, $record_file:expr) => {
        find_threshold!($c, $results, $func, $n_args:stringify!($func), $init_step, $record_file);
    };
    ($c:ident, $results:ident, $func:ident, $n_args:tt : $test_file_name:expr, $init_step:expr, $record_file:expr) => {
        'bench_block: {
            let mut group = $c.benchmark_group(format!("threshold_{}", stringify!($func)));
            let root_path = get_criterion_root();
            let func_str = stringify!($func);
            let par_path = root_path.join(format!("threshold_{}/batch", func_str));
            let ser_path = root_path.join(format!("threshold_{}/serial", func_str));
            let test_paths = &find_test_files!($func : $test_file_name);

            let mut record = Record {
                function_name: func_str.to_string(),
                threshold: $init_step,
                par_time: f64::NAN,
                ser_time: f64::NAN,
                ratio: f64::NAN,
                exit_cond: ExitCond::None,
            };

            // --- 1. Expansion phase ---
            let mut low = STEP;
            let mut high = $init_step;

            loop {
                record.threshold = ((high + STEP / 2) / STEP) * STEP;
                record.threshold = record.threshold.max(STEP);

                bench_func!(@run group, $func, test_paths, record.threshold, $n_args);

                let Ok(par_time) = extract_criterion_mean(
                    &par_path.join(record.threshold.to_string()).join("new").join("estimates.json")
                ) else {
                    group.finish();
                    break 'bench_block;
                };
                let Ok(ser_time) = extract_criterion_mean(
                    &ser_path.join(record.threshold.to_string()).join("new").join("estimates.json")
                ) else {
                    group.finish();
                    break 'bench_block;
                };
                record.par_time = par_time;
                record.ser_time = ser_time;
                record.ratio = record.par_time / record.ser_time;

                if record.ratio < 1.0 {
                    // Overshoot found
                    // --- 2. Binary search phase ---
                    let mut last_mid = record.threshold;

                    for n in 0..MAX_ITER {
                        let mid = ((low + high) / 2 / STEP) * STEP;
                        if mid == last_mid {
                            record.exit_cond = ExitCond::MinStep;
                            break;
                        }

                        bench_func!(@run group, $func, test_paths, mid, $n_args);

                        let Ok(par_time) = extract_criterion_mean(
                            &par_path.join(mid.to_string()).join("new").join("estimates.json")
                        ) else {
                            group.finish();
                            break 'bench_block;
                        };
                        let Ok(ser_time) = extract_criterion_mean(
                            &ser_path.join(mid.to_string()).join("new").join("estimates.json")
                        ) else {
                            group.finish();
                            break 'bench_block;
                        };
                        let ratio = par_time / ser_time;

                        if ratio < 1.0 - ACCEPT_THRESHOLD {
                            high = mid;
                        } else if ratio > 1.0 + ACCEPT_THRESHOLD {
                            low = mid;
                        } else {
                            record.threshold = mid;
                            record.par_time = par_time;
                            record.ser_time = ser_time;
                            record.ratio = ratio;
                            record.exit_cond = ExitCond::Converged;
                            break;
                        }

                        if high - low <= STEP {
                            record.threshold = high;
                            record.par_time = par_time;
                            record.ser_time = ser_time;
                            record.ratio = ratio;
                            record.exit_cond = ExitCond::MinStep;
                            break;
                        }

                        if n == MAX_ITER - 1 {
                            record.threshold = high;
                            record.exit_cond = ExitCond::MaxIter;
                        }

                        last_mid = mid;
                    }
                    break;
                }

                low = high;
                high *= 2;
                if high > MAX_THRESHOLD {
                    record.exit_cond = ExitCond::MaxThres;
                    break;
                }
            }
            update_or_push_record(&mut $results, record);
            save_results($record_file, &$results);
            group.finish();
        }
    };
}

macro_rules! find_thresholds {
    ($c:ident, $results:ident, $record_file:expr => { $(( $func:ident, $n_args:tt $(: $test_file_name:expr)?, $init_step:expr )),* $(,)? }) => {
        $(
            find_threshold!($c, $results, $func, $n_args $(: $test_file_name)?, $init_step, &$record_file);
        )*
    };
}

fn bench_thresholds(c: &mut Criterion) {
    let record_file = get_record_file();
    let mut results = load_results(&record_file);

    find_thresholds!(c, results, record_file => {
        (ellipk, 1, 3000),
        (ellipe, 1, 3300),
        (ellipf, 2, 400),
        (ellipeinc, 2, 300),
        (ellippi, 2, 200),
        (ellippiinc, 3, 200),
        (ellippiinc_bulirsch, 3:"ellippiinc", 300),
        (ellipd, 1, 600),
        (ellipdinc, 2, 500),
        (ellipke, 1:"ellipk", 3000),
        (cel, 4, 600),
        (cel1, 1, 1500),
        (cel2, 3, 700),
        (cel3, 2, 1400),
        (el1, 2, 600),
        (el2, 4, 600),
        (el3, 3, 500),
        (elliprf, 3, 600),
        (elliprg, 3, 500),
        (elliprj, 4, 300),
        (elliprc, 2, 800),
        (elliprd, 3, 500),
        (jacobi_zeta, 2, 200),
        (heuman_lambda, 2, 200),
    });
}

criterion_group!(benches, bench_thresholds);
criterion_main!(benches);
