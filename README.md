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

Use [Ellip-Rayon](https://github.com/p-sira/ellip/tree/main/ellip-rayon) to parallelize and improve performance for large inputs. Ellip is also available for Python via [EllipPy](https://github.com/p-sira/ellippy).

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
- Miscellaneous functions
    - `jacobi_zeta`: Jacobi Zeta function (Z). 
    - `heuman_lambda`: Heuman Lambda function (Λ0).

## Testing

In the unit tests, the functions are tested against the Boost Math and Wolfram test data. Since Ellip accepts the argument `m` (parameter) instead of `k` (modulus) to allow larger domain support, the full accuracy report uses exclusively the Wolfram data. **The full accuracy report can be found [here](https://github.com/p-sira/ellip/blob/main/tests)**, along with the test data and test generation scripts. The performance benchmark is presented to provide comparison between functions in Ellip. Comparing performance with other libraries is non-trivial, since they accept different domains.

Benchmark on AMD Ryzen 5 4600H with Radeon Graphics @3.0 GHz running x86_64-unknown-linux-gnu rustc 1.90.0 using ellip v0.5.6 at `f64` precision (ε=2.2204460492503131e-16).

### Legendre's Elliptic Integrals
| Function            | Median Error (ε) | Max Error (ε) | Mean Performance |
|---------------------|------------------|---------------|------------------|
| ellipk              | 0.00             | 108.14        | 14.8 ns          |
| ellipe              | 0.00             | 3.00          | 13.1 ns          |
| ellipf              | 0.00             | 7.47          | 98.5 ns          |
| ellipeinc           | 0.00             | 24.66         | 157.6 ns         |
| ellippi             | 0.00             | 36.35         | 166.1 ns         |
| ellippiinc          | 0.00             | 395.31        | 232.9 ns         |
| ellippiinc_bulirsch | 0.00             | 395.31        | 189.2 ns         |
| ellipd              | 0.00             | 2.64          | 30.0 ns          |
| ellipdinc           | 0.00             | 8.38          | 98.9 ns          |

### Bulirsch's Elliptic Integrals
| Function | Median Error (ε) | Max Error (ε) | Mean Performance |
|----------|------------------|---------------|------------------|
| cel      | 0.62             | 36.94         | 32.8 ns          |
| cel1     | 0.00             | 8.68          | 11.2 ns          |
| cel2     | 0.00             | 3.47          | 21.6 ns          |
| el1      | 0.00             | 1.70          | 36.7 ns          |
| el2      | 0.00             | 74.60         | 52.0 ns          |
| el3      | 0.00             | 53.21         | 103.8 ns         |

### Carlson's Symmetric Integrals
| Function | Median Error (ε) | Max Error (ε) | Mean Performance |
|----------|------------------|---------------|------------------|
| elliprf  | 0.00             | 1.57          | 46.1 ns          |
| elliprg  | 0.00             | 5.25          | 99.1 ns          |
| elliprj  | 0.56             | 136.97        | 165.5 ns         |
| elliprc  | 0.00             | 2.82          | 22.5 ns          |
| elliprd  | 0.00             | 6.25          | 75.9 ns          |

### Miscellaneous Functions
| Function      | Median Error (ε) | Max Error (ε) | Mean Performance |
|---------------|------------------|---------------|------------------|
| jacobi_zeta   | 0.00             | 8.66          | 207.7 ns         |
| heuman_lambda | 0.00             | 2.86          | 333.8 ns         |

---

Learn more at [docs.rs](https://docs.rs/ellip).