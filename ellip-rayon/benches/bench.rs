/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use std::path::{Path, PathBuf};

use criterion::{BenchmarkId, Criterion, criterion_group};
use ellip_dev_utils::{
    parser,
    test_report::{Case, format_performance},
};
use itertools::izip;
use rand::{SeedableRng, seq::IndexedRandom};
use rayon::iter::{IntoParallelRefIterator, ParallelBridge, ParallelIterator};

const SIZES: [usize; 6] = [100, 200, 300, 500, 700, 1000];
const MULTIPLE: usize = 100;

fn get_cases(test_paths: &Vec<PathBuf>, n: usize) -> Vec<Case<f64>> {
    if n % MULTIPLE != 0 {
        panic!(concat!("n must be multiple of ", stringify!(MULTIPLE), "."));
    }

    let cases = test_paths
        .iter()
        .flat_map(|test_path| parser::read_wolfram_data(test_path.to_str().unwrap()).unwrap())
        .collect::<Vec<Case<f64>>>();
    let mut rng = rand::rngs::StdRng::seed_from_u64(0);

    let base_cases: Vec<Case<f64>> = cases.choose_multiple(&mut rng, n).cloned().collect();
    let mut result = Vec::with_capacity(n);
    for _ in 0..n / MULTIPLE {
        result.extend_from_slice(&base_cases);
    }
    result
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
            $(
                let test_paths = &find_test_files!($func $($test_file_name)?);
                for &size in SIZES.iter() {
                    generate_benchmarks!(@bench group, $func, test_paths, size, $n_args);
                }
            )*
            group.finish();
        }

        macro_rules! get_func_names {
            () => {
                [$(stringify!($func)),*]
            };
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

macro_rules! get_data_paths {
    ($bench_case_dir:expr) => {
        SIZES.map(|size| {
            $bench_case_dir
                .join(size.to_string())
                .join("new")
                .join("estimates.json")
        })
    };
}

/// Criterion benchmark estimates structure
#[derive(Debug, serde::Deserialize)]
struct CriterionEstimates {
    mean: CriterionMean,
}

/// Criterion mean structure
#[derive(Debug, serde::Deserialize)]
struct CriterionMean {
    point_estimate: f64,
}

pub fn extract_criterion_mean(path: &PathBuf) -> f64 {
    use std::fs;
    let content = fs::read_to_string(path)
        .inspect_err(|_| {
            eprintln!(
                "Cannot read estimates.json file: {}",
                path.to_str().unwrap_or("nan")
            )
        })
        .unwrap();

    let estimates: CriterionEstimates = serde_json::from_str(&content)
        .inspect_err(|_| {
            eprintln!(
                "Cannot parse estimates.json file: {}",
                path.to_str().unwrap_or("nan")
            )
        })
        .unwrap();

    estimates.mean.point_estimate
}

fn print_threshold() {
    debug_assert!(SIZES.is_sorted());
    let root_path = Path::new("target/criterion/par_threshold");
    for func_name in get_func_names!() {
        let par_path = root_path.join(func_name.to_owned() + "_par");
        let ser_path = root_path.join(func_name.to_owned() + "_ser");
        let par_estimate_paths = get_data_paths!(par_path);
        let ser_estimate_paths = get_data_paths!(ser_path);
        let mut found = false;
        for (size_idx, (par_path, ser_path)) in par_estimate_paths
            .iter()
            .zip(ser_estimate_paths)
            .enumerate()
        {
            let par_estimate = extract_criterion_mean(par_path);
            let ser_estimate = extract_criterion_mean(&ser_path);
            if par_estimate < ser_estimate {
                println!(
                    "fn {func_name}: thres={} (par={} < ser={})",
                    SIZES[size_idx],
                    format_performance(&par_estimate),
                    format_performance(&ser_estimate)
                );
                found = true;
                break;
            }
        }
        if !found {
            eprintln!("fn {func_name} did not reach threshold!")
        }
    }
}

fn main() {
    // print_threshold();
    benches();
    Criterion::default().configure_from_args().final_summary();
    print_threshold();
}
