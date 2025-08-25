/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

#![cfg_attr(feature = "test_force_fail", allow(unused))]
#![allow(clippy::excessive_precision)]
//! # ELLIP
//! **Ellip** is an elliptic integral functions for Rust.
//!
//! ## Why Ellip
//! Ellip is a pure-Rust implementation of [elliptic integrals](https://dlmf.nist.gov/19).
//! This means there is no dependence on C++ libraries. Ellip also provides less common
//! functions like Bulirsch's `cel` and `el`. Some applications of elliptic integrals include
//! computing the [lengths of plane curves](https://dlmf.nist.gov/19.30), magnetic field from
//! magnets of various shapes (e.g., [cylindrical](https://doi.org/10.1016/j.jmmm.2018.02.003)),
//! [astrophysics](https://dx.doi.org/10.1088/0004-637X/696/2/1616), and [string theory](https://dx.doi.org/10.1088/1126-6708/2004/03/004).
//!
//! ## Example
//! Computing the circumference of an ellipse.
//! ```
//! use ellip::*;
//!
//! fn ellipse_length(a: f64, b: f64) -> Result<f64, StrErr> {
//!     Ok(8.0 * elliprg(0.0, a * a, b * b)?)
//! }
//!
//! let ans = ellipse_length(5.0, 3.0).unwrap();
//! ellip::util::assert_close(ans, 25.526998863398124, 1e-15);
//! ```
//!
//! # Features
//! ## Legendre's complete integrals
//! - [fn@ellipk]: Complete elliptic integral of the first kind (K).
//! - [fn@ellipe]: Complete elliptic integral of the second kind (E).
//! - [fn@ellippi]: Complete elliptic integral of the third kind (Π).
//! - [fn@ellipd]: Complete elliptic integral of Legendre's type (D).
//! ## Legendre's incomplete integrals
//! - [fn@ellipf]: Incomplete elliptic integral of the first kind (F).
//! - [fn@ellipeinc]: Incomplete elliptic integral of the second kind (E).
//! - [fn@ellippiinc]: Incomplete elliptic integral of the third kind (Π).
//! - [fn@ellippiinc_bulirsch]: Faster implementation of [fn@ellippiinc].
//! - [fn@ellipdinc]: Incomplete elliptic integral of Legendre's type (D).
//! ## Bulirsch's integrals
//! - [fn@cel]: General complete elliptic integral in Bulirsch's form.
//! - [fn@cel1]: Complete elliptic integral of the first kind in Bulirsch's form.
//! - [fn@cel2]: Complete elliptic integral of the second kind in Bulirsch's form.
//! - [fn@el1]: Incomplete elliptic integral of the first kind in Bulirsch's form.
//! - [fn@el2]: Incomplete elliptic integral of the second kind in Bulirsch's form.
//! - [fn@el3]: Incomplete elliptic integral of the third kind in Bulirsch's form.
//! ## Carlson's symmetric integrals
//! - [fn@elliprf]: Symmetric elliptic integral of the first kind (RF).
//! - [fn@elliprg]: Symmetric elliptic integral of the second kind (RG).
//! - [fn@elliprj]: Symmetric elliptic integral of the third kind (RJ).
//! - [fn@elliprc]: Degenerate elliptic integral of RF (RC).
//! - [fn@elliprd]: Degenerate elliptic integral of the third kind (RD).
//! ## Feature Flags
//! - `unstable`: Enable unstable or untested features that might be changed without notice in the future.
//! - `test_force_fail`: Used for testing only. Force tests to reach code unreachable under normal circumstances.
//!
//! # Testing
//! The function results are compared with Boost Math test data and Wolfram Engine test data.
//! The accuracy report and the test data along with the test generation scripts can
//! be found [here](https://github.com/p-sira/ellip/blob/main/tests).
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
//! - NIST Digital Library, [Chapter 19: Elliptic Integrals](https://dlmf.nist.gov/19) (Carlson, 2024).
//! - Numerical calculation of elliptic integrals and elliptic functions [I](https://link.springer.com/article/10.1007/BF01397975) (Bulirsch, 1965), [II](https://doi.org/10.1007/BF01436529) (Bulirsch, 1965), and [III](https://doi.org/10.1007/BF02165405) (Bulirsch, 1969).
//!
//! Unicode-style mathematical notation are created using [Diagon](https://github.com/ArthurSonzogni/Diagon).

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

// Utilities
mod polyeval;
use polyeval::*;
pub mod util;

#[cfg(test)]
mod test_util;
