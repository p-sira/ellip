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
    str::FromStr,
};

fn params_from_file(file_path: &str) -> Vec<Vec<f64>> {
    let file = File::open(file_path).expect("Cannot open file");
    let reader = BufReader::new(file);

    let cases: Vec<Vec<f64>> = reader
        .lines()
        .map(|line| {
            line.expect("Cannot read line")
                .split_whitespace()
                .map(|arg| f64::from_str(arg).expect(&format!("Cannot parse param {}", arg)))
                .collect()
        })
        .collect();

    cases
}

fn bench_cases<M: Measurement>(
    group: &mut BenchmarkGroup<M>,
    bench_name: &str,
    func: &dyn Fn(&[f64]) -> f64,
    cases: &Vec<Vec<f64>>,
) {
    group.bench_function(bench_name, |b| {
        b.iter(|| cases.iter().for_each(|case| assert!(!func(&case).is_nan())))
    });
}

macro_rules! generate_benchmarks {
    ($name:ident, $category:expr, $($func:ident),*) => {
        pub fn $name(c: &mut Criterion) {
            $(
                let mut group = c.benchmark_group(concat!($category, "/", stringify!($func)));
                let cases = params_from_file(concat!("tests/data/boost/", stringify!($func), "_data.txt"));

                bench_cases(&mut group, stringify!($func), &|inp| {
                    $func(inp)
                }, &cases);
                
                bench_cases(&mut group, concat!("_", stringify!($func)), &|inp| {
                    $func(inp)
                }, &cases);

                group.finish();
            )*
        }
    };
}

fn ellipk(inp: &[f64]) -> f64 {
    ellip::ellipk(inp[0]).unwrap()
}

fn _ellipk(inp: &[f64]) -> f64 {
    ellip::unchecked::_ellipk(inp[0])
}

fn ellipe(inp: &[f64]) -> f64 {
    ellip::ellipe(inp[0]).unwrap()
}

fn _ellipe(inp: &[f64]) -> f64 {
    ellip::unchecked::_ellipe(inp[0])
}

fn elliprf(inp: &[f64]) -> f64 {
    ellip::elliprf(inp[0], inp[1], inp[2]).unwrap()
}

fn _elliprf(inp: &[f64]) -> f64 {
    ellip::unchecked::_elliprf(inp[0], inp[1], inp[2])
}

generate_benchmarks!(bench_legendre, "legendre", ellipk, ellipe);
generate_benchmarks!(bench_carlson, "carlson", elliprf);

criterion_group!(benches, bench_legendre, bench_carlson);
criterion_main!(benches);
