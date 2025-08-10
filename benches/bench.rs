/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use criterion::{
    criterion_group, criterion_main, measurement::Measurement, BenchmarkGroup, Criterion,
};

mod parser {
    use ellip::StrErr;
    use ellip_dev_utils::parser::parse_wolfram_str;
    use num_traits::Float;

    /// Reads wolfram data from file and returns a vector of T
    pub fn read_wolfram_data<T: Float>(file_path: &str) -> Result<Vec<Vec<T>>, StrErr> {
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

            results.push(inputs);
        }

        Ok(results)
    }
}

fn params_from_boost_file(file_path: &str) -> Vec<Vec<f64>> {
    use std::{
        fs::File,
        io::{BufRead, BufReader},
        path::Path,
        str::FromStr,
    };

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

macro_rules! wrap_functions {
    ($($func:ident : $n_args:tt),* $(,)?) => {
        $(
            fn $func(inp: &[f64]) -> f64 {
                wrap_functions!(@call $func, inp, $n_args)
            }
        )*
    };

    (@m: $($func:ident : $n_args:tt),* $(,)?) => {
        $(
            fn $func(inp: &[f64]) -> f64 {
                wrap_functions!(@call_m $func, inp, $n_args)
            }
        )*
    };

    (@call $func:ident, $inp:ident, 1) => { ellip::$func($inp[0]).unwrap() };
    (@call $func:ident, $inp:ident, 2) => { ellip::$func($inp[0], $inp[1]).unwrap() };
    (@call $func:ident, $inp:ident, 3) => { ellip::$func($inp[0], $inp[1], $inp[2]).unwrap() };
    (@call $func:ident, $inp:ident, 4) => { ellip::$func($inp[0], $inp[1], $inp[2], $inp[3]).unwrap() };

    // Convert k to m
    (@call_m $func:ident, $inp:ident, 1) => { ellip::$func($inp[0] * $inp[0]).unwrap() };
    (@call_m $func:ident, $inp:ident, 2) => { ellip::$func($inp[0], $inp[1] * $inp[1]).unwrap() };
    (@call_m $func:ident, $inp:ident, 3) => { ellip::$func($inp[0], $inp[1], $inp[2] * $inp[2]).unwrap() };
    (@call_m $func:ident, $inp:ident, 4) => { ellip::$func($inp[0], $inp[1], $inp[2], $inp[3] * $inp[3]).unwrap() };
}

macro_rules! generate_benchmarks {
    (@params boost, $func:ident) => {
        params_from_boost_file(concat!("tests/data/boost/", stringify!($func), "_data.txt"))
    };
    (@params wolfram, $func:ident) => {
        parser::read_wolfram_data(concat!("tests/data/wolfram/", stringify!($func), "_data.csv")).unwrap()
    };
    ($bench_name:ident, [$($func:ident),*], $test_data:ident) => {
        pub fn $bench_name(c: &mut Criterion) {
            $(
                let mut group = c.benchmark_group(stringify!($bench_name));
                let cases = generate_benchmarks!(@params $test_data, $func);

                bench_cases(&mut group, stringify!($func), &|inp| {
                    $func(inp)
                }, &cases);

                group.finish();
            )*
        }
    };
    (legendre, [$($func:ident : $n_args:tt),* $(,)?], $test_data:ident $(,)?) => {
        wrap_functions! {@m: $($func : $n_args),*}
        generate_benchmarks! {legendre, [$($func),*], $test_data}
    };
    ($bench_name:ident, [$($func:ident : $n_args:tt),* $(,)?], $test_data:ident $(,)?) => {
        wrap_functions! {$($func : $n_args),*}
        generate_benchmarks! {$bench_name, [$($func),*], $test_data}
    };
}

generate_benchmarks! {legendre, [ellipk:1, ellipe:1, ellipf:2, ellipeinc:2, ellippi:2, ellippiinc:3, ellipd:1, ellipdinc: 2], boost}
generate_benchmarks! {carlson, [elliprf:3, elliprg:3, elliprj:4, elliprc:2, elliprd:3], boost}
generate_benchmarks! {bulirsch, [cel:4, cel1:1, cel2:3, el1:2, el2:4, el3:3], wolfram}

criterion_group!(benches, legendre, carlson, bulirsch);
criterion_main!(benches);
