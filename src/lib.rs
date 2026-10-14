/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

#![cfg_attr(not(any(test, feature = "std")), no_std)]
#![cfg_attr(feature = "test_force_fail", allow(unused))]
#![allow(clippy::excessive_precision)]
//! # Ellip - Elliptic integrals for Rust
//!
//! **Ellip** is a pure-Rust implementation of [elliptic integrals](https://dlmf.nist.gov/19).
//! This means there is no dependence on C++ libraries. Ellip also provides less common
//! functions like Bulirsch's `cel` and `el`. Some applications of elliptic integrals include
//! computing the [lengths of plane curves](https://dlmf.nist.gov/19.30), magnetic field from
//! magnets of various shapes (e.g., [cylindrical](https://doi.org/10.1016/j.jmmm.2018.02.003)),
//! [astrophysics](https://dx.doi.org/10.1088/0004-637X/696/2/1616), and [string theory](https://dx.doi.org/10.1088/1126-6708/2004/03/004).
//!
//! Use [Ellip-Rayon](https://github.com/p-sira/ellip/tree/main/ellip-rayon) to parallelize and improve
//! performance for large inputs. Ellip is also available for Python via [EllipPy](https://github.com/p-sira/ellippy).
//!
//! ## Example
//! Calculating the perimeter of an ellipse.
//! ```
//! use ellip::*;
//!
//! fn ellipse_perimeter(a: f64, b: f64) -> Result<f64, StrErr> {
//!     Ok(8.0 * elliprg(0.0, a * a, b * b)?)
//! }
//!
//! // Example: ellipse with semi-major axis 5, semi-minor axis 3
//! println!("{}", ellipse_perimeter(5.0, 3.0).unwrap()); // 25.526998863398124
//! ```
//!
//! # Features
//! - **[mod@legendre] — Legendre's integrals**
//!     - **Complete integrals:** the first, second, and third kinds
//!       ([fn@ellipk], [fn@ellipe], [fn@ellippi]), a combined routine ([fn@ellipke]), and
//!       Legendre's type ([fn@ellipd])
//!     - **Incomplete integrals:** the first, second, and third kinds
//!       ([fn@ellipf], [fn@ellipeinc], [fn@ellippiinc]) and Legendre's type ([fn@ellipdinc])
//! - **[mod@bulirsch] — Bulirsch's integrals**
//!     - **Complete integrals:** the general complete integral
//!       ([fn@cel]), and the first, second, and third kinds ([fn@cel1], [fn@cel2], [fn@cel3])
//!     - **Incomplete integrals:** the first, second, and third kinds
//!       ([fn@el1], [fn@el2], [fn@el3])
//! - **[mod@carlson] — Carlson's symmetric integrals:** the first, second, and third kinds
//!   ([fn@elliprf], [fn@elliprg], [fn@elliprj]), and the degenerate integrals
//!   ([fn@elliprc], [fn@elliprd])
//! - **[mod@misc] — Miscellaneous functions:** [fn@jacobi_zeta] and [fn@heuman_lambda]
//!
//! # Testing
//! The full accuracy report, test datasets, and test generation scripts are available in the
//! [tests directory](https://github.com/p-sira/ellip/tree/main/tests). The test and benchmark
//! summary is available in [TESTING.md](https://github.com/p-sira/ellip/blob/main/TESTING.md).
//! See how to reproduce the results in [REPRODUCIBILITY.md](REPRODUCIBILITY.md).
//!
//! # Acknowledgment
//! Ellip is derived from multiple mathematic libraries. We thank
//! the opensource contributors for making mathematic libraries free for all.
//! Following are the main original works used in the development of Ellip.
//! Detailed credits are available in the source code.
//! - [SciPy](https://github.com/scipy/scipy/)
//! - [Cephes Math Library](https://netlib.org/cephes/)
//! - [Boost Math Library](https://www.boost.org/doc/libs/release/libs/math/)
//! - [Russell Lab](https://github.com/cpmech/russell)
//!
//! References for original implementations are:
//! - NIST Digital Library, [Chapter 19: Elliptic Integrals](https://dlmf.nist.gov/19) (Carlson, 2025).
//! - Numerical calculation of elliptic integrals and elliptic functions [I](https://link.springer.com/article/10.1007/BF01397975) (Bulirsch, 1965), [II](https://doi.org/10.1007/BF01436529) (Bulirsch, 1965), and [III](https://doi.org/10.1007/BF02165405) (Bulirsch, 1969).
//! - NIST Digital Library, [Chapter 22 Jacobian Elliptic Functions](https://dlmf.nist.gov/22) (Reinhardt and Walker, 2025).
//!
//! Unicode-style mathematical notation are created using [Diagon](https://github.com/ArthurSonzogni/Diagon).
//!
//! # Citation
//! The paper describing Ellip is published in the [Journal of Open Source Software](https://joss.theoj.org/papers/10.21105/joss.09386). If Ellip is helpful to your work, please consider citing it:
//!
//! ```text
//! Pornsiriprasert, S., (2026). Ellip: An Elliptic Integral Library for Rust. Journal of Open Source Software, 11(118), 9386, https://doi.org/10.21105/joss.09386
//! ```
//!
//! Bibtex format:
//!
//! ```bibtex
//! @article{Pornsiriprasert2026,
//!     doi = {10.21105/joss.09386},
//!    url = {https://doi.org/10.21105/joss.09386},
//!    year = {2026},
//!    publisher = {The Open Journal},
//!    volume = {11}, number = {118}, pages = {9386},
//!    author = {Pornsiriprasert, Sira},
//!    title = {Ellip: An Elliptic Integral Library for Rust},
//!    journal = {Journal of Open Source Software} }
//! ```

num_lazy::declare_nums! {@constant T}
num_lazy::declare_nums! {@special T}

mod crate_util;

/// Static error str
pub type StrErr = &'static str;

pub mod legendre;
// Legendre's complete integrals
pub use legendre::ellipd;
pub use legendre::ellipe;
pub use legendre::ellipk;
pub use legendre::ellipke;
pub use legendre::ellippi;

// Legendre's incomplete integrals
pub use legendre::ellipdinc;
pub use legendre::ellipeinc;
pub use legendre::ellipf;
pub use legendre::ellippiinc;
pub use legendre::ellippiinc_bulirsch;

// Bulirsch's integrals
pub mod bulirsch;
pub use bulirsch::cel;
pub use bulirsch::cel1;
pub use bulirsch::cel2;
pub use bulirsch::cel3;
pub use bulirsch::el1;
pub use bulirsch::el2;
pub use bulirsch::el3;

// Carlson's symmetric integrals
pub mod carlson;
pub use carlson::elliprc;
pub use carlson::elliprd;
pub use carlson::elliprf;
pub use carlson::elliprg;
pub use carlson::elliprj;

// Miscellaneous functions
pub mod misc;
pub use misc::heuman_lambda;
pub use misc::jacobi_zeta;

// Utilities
mod polyeval;
use polyeval::polyeval;
pub mod util;

#[cfg(test)]
mod test_util;
