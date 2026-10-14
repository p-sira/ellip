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

## Quick Start

Start by installing Ellip, which is `no_std` by default.

```shell
>> cargo add ellip
```

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
- **Legendre's complete integrals:** the first, second, and third kinds (`ellipk`, `ellipe`, `ellippi`), a combined routine (`ellipke`), and Legendre's type (`ellipd`)
- **Legendre's incomplete integrals:** the first, second, and third kinds (`ellipf`, `ellipeinc`, `ellippiinc`) and Legendre's type (`ellipdinc`)
- **Bulirsch's complete integrals:** the general complete integral (`cel`), and the first, second, and third kinds (`cel1`, `cel2`, `cel3`)
- **Bulirsch's incomplete integrals:** the first, second, and third kinds (`el1`, `el2`, `el3`)
- **Carlson's symmetric integrals:** the first, second, and third kinds (`elliprf`, `elliprg`, `elliprj`), and the degenerate integrals (`elliprc`, `elliprd`) 
- **Miscellaneous functions:** `jacobi_zeta` and `heuman_lambda`

## Testing

The full accuracy report, test datasets, and test generation scripts are available in the [tests directory](https://github.com/p-sira/ellip/tree/main/tests). The test and benchmark summary is available in [TESTING.md](https://github.com/p-sira/ellip/blob/main/TESTING.md). See how to reproduce the results in [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

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