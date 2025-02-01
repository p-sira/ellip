/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

#![allow(clippy::excessive_precision)]

//! # ELLIP
//! **Ellip** is an elliptic integral functions for Rust.
//!
//! # Features
//! ## Legendre's complete integrals
//! - [fn@ellipk]: Complete elliptic integral of the first kind.
//! - [fn@ellipe]: Complete elliptic integral of the second kind.
//! ## Legendre's incomplete integrals
//! - [fn@ellipf]: Incomplete elliptic integral of the first kind.
//! - [fn@ellipeinc]: Incomplete elliptic integral of the second kind.
//! ## Bulirsch's integrals
//! - [fn@cel]: General complete elliptic integral
//! ## Carlson's symmetric integrals
//! - [fn@elliprf]: Symmetric elliptic integral of the first kind.
//! - [fn@elliprg]: Symmetric elliptic integral of the second kind.
//! - [fn@elliprj]: Symmetric elliptic integral of the third kind.
//! - [fn@elliprc]: Degenerate elliptic integral of RF.
//! - [fn@elliprd]: Degenerate elliptic integral of the third kind.
//!
//! # Acknowledgment
//! Ellip is derived from multiple mathematic libraries. We thank
//! the opensource contributors for making mathematic libraries free for all.
//! Following are the main original works used in the development of Ellip.
//! Detailed credits are available in the source code.
//! - [Scipy](https://github.com/scipy/scipy/)
//! - [Cephes Math Library](https://netlib.org/cephes/)
//! - [Boost Math Library](https://www.boost.org/doc/libs/release/libs/math/)
//! - [Russell Lab](https://github.com/cpmech/russell)
//!
//! Primary mathematical reference is [Chapter 19](https://dlmf.nist.gov/19) of the NIST Digital Library
//! of Mathematical Functions, authored by [Carlson](https://dlmf.nist.gov/about/bio/BCCarlson), the legendary
//! mathematician who discovered the symmetric integrals!
//!
//! Unicode-style mathematical notation are created using an awesome tool called
//! [Diagon](https://github.com/ArthurSonzogni/Diagon).

// Legendre's complete integrals
mod ellipk;
pub use ellipk::ellipk;
mod ellipe;
pub use ellipe::ellipe;

// Legendre's incomplete integrals
mod ellipf;
pub use ellipf::ellipf;
mod ellipeinc;
pub use ellipeinc::ellipeinc;

// Bulirsch's integrals
mod cel;
pub use cel::cel;

// Carlson's symmetric integrals
mod elliprf;
pub use elliprf::elliprf;
mod elliprg;
pub use elliprg::elliprg;
mod elliprj;
pub use elliprj::elliprj;
mod elliprc;
pub use elliprc::elliprc;
mod elliprd;
pub use elliprd::elliprd;

// Utilities
mod polyeval;
use polyeval::*;
const TINY: f64 = 5.0 * f64::MIN_POSITIVE;
const BIG: f64 = 0.2 * f64::MAX;

#[cfg(test)]
mod test_util;
