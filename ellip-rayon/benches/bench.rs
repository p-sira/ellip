/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use std::io::Write;
use std::path::{Path, PathBuf};

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use ellip_dev_utils::{
    benchmark::extract_criterion_mean,
    parser,
    test_report::{Case, format_float, format_performance},
};
use itertools::izip;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};
use tabled::{Table, Tabled, settings::Style};

const MAX_THRESHOLD: usize = 10000;
const ACCEPT_THRESHOLD: f64 = 0.1;
const STEP: usize = 100;
const INIT_STEP: usize = 10 * STEP;

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

fn get_cases(test_paths: &Vec<PathBuf>, n: usize) -> Vec<Case<f64>> {
    if n % STEP != 0 {
        panic!(concat!("n must be multiple of ", stringify!(MULTIPLE), "."));
    }

    let cases = test_paths
        .iter()
        .flat_map(|test_path| parser::read_wolfram_data(test_path.to_str().unwrap()).unwrap())
        .collect::<Vec<Case<f64>>>();
    let base_cases: Vec<Case<f64>> = linspace(0, cases.len(), FILE_SAMPLE_SIZE)
        .iter()
        .map(|&i| cases[i].clone())
        .collect();
    if n > FILE_SAMPLE_SIZE {
        let mut result = Vec::with_capacity(n);
        for _ in 0..n / STEP {
            result.extend_from_slice(&base_cases);
        }
        result
    } else {
        base_cases.iter().take(n).cloned().collect()
    }
}

fn extract_params(cases: &[Case<f64>], param_index: usize) -> Vec<f64> {
    cases.iter().map(|case| case.inputs[param_index]).collect()
}

pub fn find_test_files(function_name: &str, dir: &str) -> Vec<PathBuf> {
    // Also try the wolfram directory specifically as fallback
    let pattern = format!("../tests/data/{dir}/{function_name}_*.csv");
    glob::glob(&pattern)
        .unwrap()
        .map(|path| path.unwrap())
        .collect()
}

macro_rules! find_test_files {
    ($func:ident) => {
        find_test_files(stringify!($func), "wolfram")
    };
    ($func:ident $test_file_name:expr) => {
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

macro_rules! generate_benchmarks {
    (@bench_inner $group:expr, $func:ident, $size:expr, $ser:expr, $par:expr $(,)?) => {
            $group.bench_with_input(
                BenchmarkId::new(concat![stringify!($func), "_ser"], $size),
                &$size,
                |b, _| {
                    b.iter(|| {
                        let ans: Vec<f64> = $ser;
                        assert!(!ans.is_empty())
                    });
                },
            );

            $group.bench_with_input(
                BenchmarkId::new(concat![stringify!($func), "_par"], $size),
                &$size,
                |b, _| {
                    b.iter(|| {
                        let ans: Vec<f64> = $par;
                        assert!(!ans.is_empty())
                    });
                },
            );
    };
    (@bench $group:expr, $func:ident, $test_paths:expr, $size:expr, 1) => {{
        let cases = get_cases($test_paths, $size);
        let a = extract_params(&cases, 0);
        generate_benchmarks!(
            @bench_inner $group, $func, $size,
            a.iter().map(|&x| ellip::$func(x).unwrap()).collect(),
            a.par_iter().map(|&x| ellip::$func(x).unwrap()).collect(),
        );
    }};
    (@bench $group:expr, $func:ident, $test_paths:expr, $size:expr, 2) => {{
        let cases = get_cases($test_paths, $size);
        let a = extract_params(&cases, 0);
        let b = extract_params(&cases, 1);
        generate_benchmarks!(
            @bench_inner $group, $func, $size,
            izip!(&a, &b).map(|(&x, &y)| ellip::$func(x, y).unwrap()).collect(),
            par_zip!(&a, &b).map(|(&x, &y)| ellip::$func(x, y).unwrap()).collect(),
        );
    }};
    (@bench $group:expr, $func:ident, $test_paths:expr, $size:expr, 3) => {{
        let cases = get_cases($test_paths, $size);
        let a = extract_params(&cases, 0);
        let b = extract_params(&cases, 1);
        let c = extract_params(&cases, 2);
        generate_benchmarks!(
            @bench_inner $group, $func, $size,
            izip!(&a, &b, &c).map(|(&x, &y, &z)| ellip::$func(x, y, z).unwrap()).collect(),
            par_zip!(&a, &b, &c).map(|(&x, &y, &z)| ellip::$func(x, y, z).unwrap()).collect(),
        );
    }};
    (@bench $group:expr, $func:ident, $test_paths:expr, $size:expr, 4) => {{
        let cases = get_cases($test_paths, $size);
        let a = extract_params(&cases, 0);
        let b = extract_params(&cases, 1);
        let c = extract_params(&cases, 2);
        let d = extract_params(&cases, 3);
        generate_benchmarks!(
            @bench_inner $group, $func, $size,
            izip!(&a, &b, &c, &d).map(|(&x, &y, &z, &p)| ellip::$func(x, y, z, p).unwrap()).collect(),
            par_zip!(&a, &b, &c, &d).map(|(&x, &y, &z, &p)| ellip::$func(x, y, z, p).unwrap()).collect(),
        );
    }};
    ($($func:ident : $n_args:tt $(: $test_file_name:expr)?),* $(,)?) => {
        pub fn par_threshold(c: &mut Criterion) {
            let mut group = c.benchmark_group("par_threshold");
            let root_path = Path::new("target/criterion/par_threshold");
            let record_file = Path::new("benches/par_threshold.md");
            let mut record_file = std::fs::File::create(record_file).unwrap();
            let mut results = Vec::new();
            $({
                let test_paths = &find_test_files!($func $($test_file_name)?);
                let func_name = stringify!($func);
                let par_path = root_path.join(func_name.to_owned() + "_par");
                let ser_path = root_path.join(func_name.to_owned() + "_ser");

                let mut record = Record {
                    function_name: func_name.to_string(),
                    threshold: INIT_STEP,
                    par_time: f64::NAN,
                    ser_time: f64::NAN,
                    ratio: f64::NAN,
                    exit_cond: ExitCond::None,
                };

                // --- 1. Expansion phase ---
                let mut low = STEP;
                let mut high = INIT_STEP;

                loop {
                    // Align to STEP
                    record.threshold = ((high + STEP / 2) / STEP) * STEP;
                    record.threshold = record.threshold.max(STEP);

                    generate_benchmarks!(@bench group, $func, test_paths, record.threshold, $n_args);
                    record.par_time = extract_criterion_mean(
                        &par_path.join(record.threshold.to_string()).join("new").join("estimates.json")
                    ).unwrap();
                    record.ser_time = extract_criterion_mean(
                        &ser_path.join(record.threshold.to_string()).join("new").join("estimates.json")
                    ).unwrap();
                    record.ratio = record.par_time / record.ser_time;

                    if record.ratio < 1.0 {
                        // Overshoot found
                        // --- 2. Binary search phase ---
                        for n in 0..MAX_ITER {
                            let mid = ((low + high) / 2 / STEP) * STEP;
                            generate_benchmarks!(@bench group, $func, test_paths, mid, $n_args);
                            let par_time = extract_criterion_mean(
                                &par_path.join(mid.to_string()).join("new").join("estimates.json")
                            ).unwrap();
                            let ser_time = extract_criterion_mean(
                                &ser_path.join(mid.to_string()).join("new").join("estimates.json")
                            ).unwrap();
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

                results.push(record);
            })*
            group.finish();
            let result_str = Table::new(results).with(Style::markdown()).to_string();
            record_file.write_all(result_str.as_bytes()).unwrap();
        }
    };
}

generate_benchmarks!(
    ellipk:1, ellipe:1, ellippi:2, ellipd:1,
    ellipf:2, ellipeinc:2, ellippiinc:3, ellippiinc_bulirsch:3:"ellippiinc", ellipdinc:2,
    elliprf:3, elliprg:3, elliprj:4, elliprc:2, elliprd:3,
    cel:4, cel1:1, cel2:3, el1:2, el2:4, el3:3,
    jacobi_zeta:2, heuman_lambda:2,
);

criterion_group!(benches, par_threshold);
criterion_main!(benches);
