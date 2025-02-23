/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use criterion::{
    criterion_group, criterion_main, measurement::Measurement, BenchmarkGroup, Criterion,
};
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
    str::FromStr,
};

fn params_from_file(file_path: &str) -> Vec<Vec<f64>> {
    let path = Path::new(file_path);
    if !path.exists() {
        panic!("Test data not found: {}", file_path);
    }

    let file = File::open(file_path).expect("Cannot open file");
    let reader = BufReader::new(file);

    let cases: Vec<Vec<f64>> = reader
        .lines()
        .map(|line| {
            line.expect("Cannot read line")
                .split_whitespace()
                .map(|arg| {
                    f64::from_str(arg).unwrap_or_else(|_| panic!("Cannot parse param {}", arg))
                })
                .collect()
        })
        .collect();

    cases
}

fn bench_cases<M: Measurement>(
    group: &mut BenchmarkGroup<M>,
    bench_name: &str,
    func: &dyn Fn(&[f64]) -> f64,
    cases: &[Vec<f64>],
) {
    group.bench_function(bench_name, |b| {
        b.iter(|| cases.iter().for_each(|case| assert!(!func(case).is_nan())))
    });
}

macro_rules! generate_benchmarks {
    ($name:ident, $category:expr, $($func:ident),*) => {
        pub fn $name(c: &mut Criterion) {
            $(
                let mut group = c.benchmark_group($category);
                let cases = params_from_file(concat!("tests/data/boost/", stringify!($func), "_data.txt"));

                bench_cases(&mut group, stringify!($func), &|inp| {
                    $func(inp)
                }, &cases);

                group.finish();
            )*
        }
    };
}

fn ellipk(inp: &[f64]) -> f64 {
    ellip::ellipk(inp[0] * inp[0]).unwrap()
}

fn ellipe(inp: &[f64]) -> f64 {
    ellip::ellipe(inp[0] * inp[0]).unwrap()
}

fn ellipf(inp: &[f64]) -> f64 {
    ellip::ellipf(inp[0], inp[1] * inp[1]).unwrap()
}

fn ellipeinc(inp: &[f64]) -> f64 {
    ellip::ellipeinc(inp[0], inp[1] * inp[1]).unwrap()
}

fn ellippi(inp: &[f64]) -> f64 {
    ellip::ellippi(inp[0], inp[1] * inp[1]).unwrap()
}

fn ellippiinc(inp: &[f64]) -> f64 {
    ellip::ellippiinc(inp[1], inp[0], inp[2] * inp[2]).unwrap()
}

fn elliprf(inp: &[f64]) -> f64 {
    ellip::elliprf(inp[0], inp[1], inp[2]).unwrap()
}

fn elliprg(inp: &[f64]) -> f64 {
    ellip::elliprg(inp[0], inp[1], inp[2]).unwrap()
}

fn elliprj(inp: &[f64]) -> f64 {
    ellip::elliprj(inp[0], inp[1], inp[2], inp[3]).unwrap()
}

fn elliprc(inp: &[f64]) -> f64 {
    ellip::elliprc(inp[0], inp[1]).unwrap()
}

fn elliprd(inp: &[f64]) -> f64 {
    ellip::elliprd(inp[0], inp[1], inp[2]).unwrap()
}

generate_benchmarks!(
    bench_legendre,
    "legendre",
    ellipk,
    ellipe,
    ellipf,
    ellipeinc,
    ellippi,
    ellippiinc
);

generate_benchmarks!(
    bench_carlson,
    "carlson",
    elliprf,
    elliprg,
    elliprj,
    elliprc,
    elliprd
);

criterion_group!(benches, bench_legendre, bench_carlson);
criterion_main!(benches);
