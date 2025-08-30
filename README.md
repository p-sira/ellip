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
        <img src="https://codecov.io/github/p-sira/ellip/graph/badge.svg?token=JVM89PIP5K"/> 
    </a>
</p>

<big><p align="center"> 
Elliptic integrals for Rust 
</p></big>

Ellip is a pure-Rust implementation of [elliptic integrals](https://dlmf.nist.gov/19). Ellip also provides less common functions like Bulirsch's `cel` and `el`. Some applications of the elliptic integrals include computing the [lengths of plane curves](https://dlmf.nist.gov/19.30), [magnetism](https://doi.org/10.1016/j.jmmm.2018.02.003), [astrophysics](https://dx.doi.org/10.1088/0004-637X/696/2/1616), and [string theory](https://dx.doi.org/10.1088/1126-6708/2004/03/004).

## Quick Start

Start by installing Ellip.
```shell
>> cargo add ellip
```

Let's compute the circumference of an ellipse.

```rust
use ellip::*;

fn ellipse_length(a: f64, b: f64) -> Result<f64, StrErr> {
    Ok(8.0 * elliprg(0.0, a * a, b * b)?)
}

let ans = ellipse_length(5.0, 3.0).unwrap();
ellip::util::assert_close(ans, 25.526998863398124, 1e-15);
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

# Testing

In the unit tests, the functions are tested against the Boost Math and Wolfram test data. Since Ellip accepts the argument `m` (parameter) instead of `k` (modulus) to allow larger domain support, the full accuracy report uses exclusively the Wolfram data. The full accuracy report, test data, and test generation scripts can be found [here](https://github.com/p-sira/ellip/blob/main/tests). The performance benchmark is presented to provide comparison between functions in Ellip. Comparing performance with other libraries is non-trivial, since they accept different domains.

Benchmark on AMD Ryzen 5 4600H with Radeon Graphics running x86_64-unknown-linux-gnu rustc 1.89.0 using ellip v0.4.0 at `f64` precision (ε≈2.22e-16).

### Legendre's Elliptic Integrals
| Function            | Median Error (ε) | Max Error (ε) | Mean Performance |
|---------------------|------------------|---------------|------------------|
| ellipk              | 0.00             | 108.14        | 14.8 ns          |
| ellipe              | 0.00             | 3.00          | 13.3 ns          |
| ellipf              | 0.66             | 7.47          | 104.1 ns         |
| ellipeinc           | 0.70             | 24.66         | 168.1 ns         |
| ellippi             | 0.53             | 36.35         | 170.1 ns         |
| ellippiinc          | 0.78             | 1.04e3        | 279.6 ns         |
| ellippiinc_bulirsch | 0.83             | 1.04e3        | 227.7 ns         |
| ellipd              | 0.60             | 2.64          | 30.5 ns          |
| ellipdinc           | 1.00             | 8.38          | 106.1 ns         |

### Bulirsch's Elliptic Integrals
| Function | Median Error (ε) | Max Error (ε) | Mean Performance |
|----------|------------------|---------------|------------------|
| cel      | 0.70             | 38.34         | 34.7 ns          |
| cel1     | 0.00             | 8.68          | 11.7 ns          |
| cel2     | 0.61             | 3.97          | 23.1 ns          |
| el1      | 0.00             | 1.60          | 38.4 ns          |
| el2      | 0.70             | 79.92         | 57.3 ns          |
| el3      | 0.70             | 46.32         | 120.6 ns         |

### Carlson's Symmetric Integrals
| Function | Median Error (ε) | Max Error (ε) | Mean Performance |
|----------|------------------|---------------|------------------|
| elliprf  | 0.00             | 1.75          | 47.3 ns          |
| elliprg  | 0.00             | 2.45          | 100.1 ns         |
| elliprj  | 0.67             | 5.42e7        | 166.5 ns         |
| elliprc  | 0.00             | 2.82          | 23.2 ns          |
| elliprd  | 0.62             | 6.49          | 76.1 ns          |

---

Learn more at [docs.rs](https://docs.rs/ellip).