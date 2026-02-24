<h1 align="center">
    <a href="https://github.com/p-sira/ellip/">
        <img src="https://github.com/p-sira/ellip/blob/main/logo/ellip-logo.svg?raw=true" alt="ELLIP" width="300">
    </a>
</h1>

<p align="center">
    <a href="https://opensource.org/license/BSD-3-clause">
        <img src="https://img.shields.io/badge/License-BSD--3--Clause-brightgreen.svg" alt="License">
    </a>
    <a href="https://crates.io/crates/ellip">
        <img src="https://img.shields.io/crates/v/ellip" alt="Crate">
    </a>
    <a href="https://crates.io/crates/ellip">
        <img src="https://img.shields.io/crates/d/ellip" alt="Total Downloads">
    </a>
    <a href="https://docs.rs/ellip">
        <img src="https://img.shields.io/badge/Docs-docs.rs-blue" alt="Documentation">
    </a>
    <a href="https://codecov.io/github/p-sira/ellip" > 
        <img src="https://codecov.io/github/p-sira/ellip/graph/badge.svg?token=JVM89PIP5K"> 
    </a>
    <a style="border-width:0" href="https://doi.org/10.21105/joss.09386">
        <img src="https://joss.theoj.org/papers/10.21105/joss.09386/status.svg" alt="DOI badge">
    </a>
</p>

<big><p align="center"> 
Elliptic integrals for Rust 
</p></big>

Ellip is a pure-Rust implementation of [elliptic integrals](https://dlmf.nist.gov/19). Ellip also provides less common functions like Bulirsch's `cel` and `el`. Some applications of the elliptic integrals include computing the [lengths of plane curves](https://dlmf.nist.gov/19.30), [magnetism](https://doi.org/10.1016/j.jmmm.2018.02.003), [astrophysics](https://dx.doi.org/10.1088/0004-637X/696/2/1616), and [string theory](https://dx.doi.org/10.1088/1126-6708/2004/03/004).

Use [Ellip-Rayon](https://github.com/p-sira/ellip/tree/main/ellip-rayon) to parallelize and improve performance for large inputs. Ellip is also available for Python via [EllipPy](https://github.com/p-sira/ellippy).

## Installation

Start by installing Ellip.

```shell
>> cargo add ellip
```

By default, Ellip installs with `no_std` and uses `libm` as the mathematical backend.
To use `std`, install Ellip using the following command:

```shell
>> cargo add ellip --no-default-features --features=std
```

## Quick Start

Let's compute the perimeter of an ellipse.

```rust
use ellip::*;

fn ellipse_perimeter(a: f64, b: f64) -> Result<f64, StrErr> {
    Ok(8.0 * elliprg(0.0, a * a, b * b)?)
}

// Example: ellipse with semi-major axis 5, semi-minor axis 3
println!("{}", ellipse_perimeter(5.0, 3.0).unwrap()); // 25.526998863398124
```

Learn more at [doc.rs](https://docs.rs/ellip).

## Features
- Legendre's complete integrals
    - `ellipk`: Complete elliptic integral of the first kind (K).
    - `ellipe`: Complete elliptic integral of the second kind (E).
    - `ellippi`: Complete elliptic integral of the third kind (Π).
    - `ellipd`: Complete elliptic integral of Legendre's type (D).
- Legendre's incomplete integrals
    - `ellipf`: Incomplete elliptic integral of the first kind (F).
    - `ellipeinc`: Incomplete elliptic integral of the second kind (E).
    - `ellippiinc`: Incomplete elliptic integral of the third kind (Π).
    - `ellipdinc`: Incomplete elliptic integral of Legendre's type (D).
- Bulirsch's integrals
    - `cel`: General complete elliptic integral in Bulirsch's form.
    - `cel1`: Complete elliptic integral of the first kind in Bulirsch's form.
    - `cel2`: Complete elliptic integral of the second kind in Bulirsch's form.
    - `el1`: Incomplete elliptic integral of the first kind in Bulirsch's form.
    - `el2`: Incomplete elliptic integral of the second kind in Bulirsch's form.
    - `el3`: Incomplete elliptic integral of the third kind in Bulirsch's form.
- Carlson's symmetric integrals
    - `elliprf`: Symmetric elliptic integral of the first kind (RF).
    - `elliprg`: Symmetric elliptic integral of the second kind (RG).
    - `elliprj`: Symmetric elliptic integral of the third kind (RJ).
    - `elliprc`: Degenerate elliptic integral of RF (RC).
    - `elliprd`: Degenerate elliptic integral of the third kind (RD).
- Miscellaneous functions
    - `jacobi_zeta`: Jacobi Zeta function (Z). 
    - `heuman_lambda`: Heuman Lambda function (Λ0).

## Testing

In the unit tests, the functions are tested against the Boost Math and Wolfram test data. Since Ellip accepts the argument `m` (parameter) instead of `k` (modulus) to allow larger domain support, the full accuracy report uses exclusively the Wolfram data. **The full accuracy report can be found [here](https://github.com/p-sira/ellip/blob/main/tests)**, along with the test data and test generation scripts. The performance benchmark is presented to provide comparison between functions in Ellip. Comparing performance with other libraries is non-trivial, since they accept different domains of input.

{{TESTSUMMARY}}

## Reproducibility

This section describes how to reproduce the accuracy reports, test datasets, benchmarks, figures, and tables.

### Setup the project

First, clone the repository:

```sh
git clone https://github.com/p-sira/ellip.git
cd ellip
```

Then, build the project:

```sh
cargo build --workspace
```

### Run tests

For detailed information on tests, see [tests/README.md](https://github.com/p-sira/ellip/blob/main/tests/README.md).

### Benchmark

Ellip's benchmark collects the test files associated with each function and reports the total execution time:

```sh
cargo bench
```

This produces raw benchmark output under `target/criterion/`. Note that the results shown in the README are normalized to per-function call averages. See the [Generate tables](#generate-tables) section for details on generating the summary tables.

### Generate tables

To generate the accuracy table in the [test report](https://github.com/p-sira/ellip/blob/main/tests):

```sh
cargo run --example generate_error_report
```

This compares Ellip's results against Wolfram test data and generates the accuracy report.

To generate the test and benchmark summary table as shown in the README, first run `cargo bench` to collect benchmark data. Then run:

```sh
cargo run --example generate_test_summary
```

This script compares results against Wolfram data, extracts benchmark results from `target/criterion/`, normalizes them to average time per function call, and summarizes everything in a single table.

### Generate figures

To generate function plots:

```sh
cargo run -p ellip-plot-graph --bin [function-name]
```

See available plots in [ellip-plot-graph/src/bin](https://github.com/p-sira/ellip/blob/main/ellip-plot-graph/src/bin/)

## Citation

The paper describing Ellip is published in the [Journal of Open Source Software](https://joss.theoj.org/papers/10.21105/joss.09386). If Ellip is helpful to your work, please consider citing it:

```text
Pornsiriprasert, S., (2026). Ellip: An Elliptic Integral Library for Rust. Journal of Open Source Software, 11(118), 9386, https://doi.org/10.21105/joss.09386
```

Bibtex format:

```bibtex
@article{Pornsiriprasert2026,
    doi = {10.21105/joss.09386}, 
    url = {https://doi.org/10.21105/joss.09386}, 
    year = {2026}, 
    publisher = {The Open Journal}, 
    volume = {11}, number = {118}, pages = {9386},
    author = {Pornsiriprasert, Sira}, 
    title = {Ellip: An Elliptic Integral Library for Rust}, 
    journal = {Journal of Open Source Software} } 
```

---

Learn more at [docs.rs](https://docs.rs/ellip).