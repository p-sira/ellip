# Reproducibility

This section describes how to reproduce the accuracy reports, test datasets, benchmarks, figures, and tables.

## Setup the project

First, clone the repository:

```sh
git clone https://github.com/p-sira/ellip.git
cd ellip
```

Then, build the project:

```sh
cargo build --workspace
```

## Run tests

For detailed information on tests, see [tests/README.md](tests/README.md).

## Benchmark

Ellip's benchmark collects the test files associated with each function and reports the total execution time:

```sh
cargo bench
```

This produces raw benchmark output under `target/criterion/`. Note that the results shown in [TESTING.md](TESTING.md) are normalized to per-function call averages. See the [Generate tables](#generate-tables) section for details on generating the summary tables.

## Generate tables

To generate the accuracy table in the [test report](tests):

```sh
cargo run --example generate_error_report
```

This compares Ellip's results against Wolfram test data and generates the accuracy report.

To generate the test and benchmark summary table, first run `cargo bench` to collect benchmark data. Then run:

```sh
cargo run --example generate_test_summary
```

This script compares results against Wolfram data, extracts benchmark results from `target/criterion/`, normalizes them to average time per function call, and summarizes everything in a single table.

## Generate figures

To generate function plots:

```sh
cargo run -p ellip-plot-graph --bin [function-name]
```

See available plots in [ellip-plot-graph/src/bin](ellip-plot-graph/src/bin/)
