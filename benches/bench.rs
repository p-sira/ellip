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
    (@find_test $func:ident) => {
        ellip_dev_utils::file::find_test_files(stringify!($func), "wolfram")
    };
    (@find_test $func:ident $test_file_name:expr) => {
        ellip_dev_utils::file::find_test_files($test_file_name, "wolfram")
    };
    ($bench_name:ident, [$($func:ident : $n_args:tt $(: $test_file_name:expr)?),* $(,)?] $(,)?) => {
        wrap_functions! {$($func : $n_args),*}
        pub fn $bench_name(c: &mut Criterion) {
            $(
                let mut group = c.benchmark_group(stringify!($bench_name));

                let test_paths = generate_benchmarks!(@find_test $func $($test_file_name)?);
                if test_paths.is_empty() {
                    eprintln!("No test files found for function: {}", stringify!($func));
                    return;
                }

                let cases = test_paths
                    .iter()
                    .flat_map(|test_path| parser::read_wolfram_data(test_path.to_str().unwrap()).unwrap())
                    .collect::<Vec<Vec<f64>>>();

                bench_cases(&mut group, stringify!($func), &|inp| {
                    $func(inp)
                }, &cases);

                group.finish();
            )*
        }
    };
}

generate_benchmarks! {legendre, [ellipk:1, ellipe:1, ellipf:2, ellipeinc:2, ellippi:2, ellippiinc:3, ellippiinc_bulirsch:3:"ellippiinc", ellipd:1, ellipdinc: 2]}
generate_benchmarks! {carlson, [elliprf:3, elliprg:3, elliprj:4, elliprc:2, elliprd:3]}
generate_benchmarks! {bulirsch, [cel:4, cel1:1, cel2:3, el1:2, el2:4, el3:3]}
generate_benchmarks! {misc, [jacobi_zeta:2, heuman_lambda:2]}

criterion_group!(benches, legendre, carlson, bulirsch, misc);
criterion_main!(benches);
