/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use std::path::{Path, PathBuf};

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use ellip_dev_utils::{
    benchmark::extract_criterion_mean,
    parser,
    test_report::{format_performance, Case},
};
use itertools::izip;
use rayon::iter::{IntoParallelRefIterator, ParallelBridge, ParallelIterator};

const ACCEPT_THRESHOLD: f64 = 0.1;
const INITIAL_STEP_SIZE: f64 = 2.0; // Initial step size as fraction
const DECAY_FACTOR: f64 = 0.8; // Exponential decay factor
const MIN_STEP_SIZE: f64 = 0.01; // Minimum step size as fraction
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
            izip!(&a, &b).par_bridge().map(|(&x, &y)| ellip::$func(x, y).unwrap()).collect(),
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
            izip!(&a, &b, &c).par_bridge().map(|(&x, &y, &z)| ellip::$func(x, y, z).unwrap()).collect(),
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
            izip!(&a, &b, &c, &d).par_bridge().map(|(&x, &y, &z, &p)| ellip::$func(x, y, z, p).unwrap()).collect(),
        );
    }};
    ($($func:ident : $n_args:tt $(: $test_file_name:expr)?),* $(,)?) => {
        pub fn par_threshold(c: &mut Criterion) {
            let mut group = c.benchmark_group("par_threshold");
            let root_path = Path::new("target/criterion/par_threshold");
            $({
                let test_paths = &find_test_files!($func $($test_file_name)?);
                let func_name = stringify!($func);
                let par_path = root_path.join(func_name.to_owned() + "_par");
                let ser_path = root_path.join(func_name.to_owned() + "_ser");

                let mut size = STEP;
                let mut step_size = INITIAL_STEP_SIZE;
                let mut visited: std::collections::HashSet<usize> = std::collections::HashSet::new();
                let mut last_direction = 0;

                for iter in 0..MAX_ITER {
                    // Round size to nearest STEP multiple
                    size = ((size + STEP / 2) / STEP) * STEP;
                    size = size.max(STEP); // Ensure minimum size

                    if visited.contains(&size) {
                        // We've been here before, reduce step size
                        step_size *= DECAY_FACTOR;
                        if step_size < MIN_STEP_SIZE {
                            println!(
                                "fn {func_name}: thres={} (converged after {} iterations)",
                                size, iter + 1
                            );
                            break;
                        }
                        continue;
                    }

                    visited.insert(size);

                    generate_benchmarks!(@bench group, $func, test_paths, size, $n_args);

                    let par_estimate = extract_criterion_mean(&par_path.join(size.to_string()).join("new").join("estimates.json")).unwrap();
                    let ser_estimate = extract_criterion_mean(&ser_path.join(size.to_string()).join("new").join("estimates.json")).unwrap();
                    let ratio = par_estimate / ser_estimate;

                    if ratio < 1.0 - ACCEPT_THRESHOLD {
                        // Parallel is significantly faster, try smaller size
                        let current_direction = -1;
                        if last_direction == 1 {
                            // Direction changed, apply decay
                            step_size *= DECAY_FACTOR;
                        }
                        last_direction = current_direction;

                        let step_amount = (size as f64 * step_size) as usize;
                        let step_amount = step_amount.max(STEP); // Minimum step
                        size = size.saturating_sub(step_amount);

                    } else if ratio > 1.0 + ACCEPT_THRESHOLD {
                        // Serial is significantly faster, try larger size
                        let current_direction = 1;
                        if last_direction == -1 {
                            // Direction changed, apply decay
                            step_size *= DECAY_FACTOR;
                        }
                        last_direction = current_direction;

                        let step_amount = (size as f64 * step_size) as usize;
                        let step_amount = step_amount.max(STEP); // Minimum step
                        size += step_amount;

                    } else {
                        // Within acceptable threshold
                        println!(
                            "fn {func_name}: thres={} (par={} ≈ ser={}, ratio={:.3})",
                            size,
                            format_performance(&par_estimate),
                            format_performance(&ser_estimate),
                            ratio
                        );
                        break;
                    }

                    // Apply minimum step size check
                    if step_size < MIN_STEP_SIZE {
                        println!(
                            "fn {func_name}: thres={} (step size too small, ratio={:.3})",
                            size, ratio
                        );
                        break;
                    }

                    if iter == MAX_ITER - 1 {
                        eprintln!("fn {func_name}: failed to find threshold within max iterations, last ratio={:.3}", ratio);
                    }
                }
            })*
            group.finish();
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
