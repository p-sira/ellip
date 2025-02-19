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
//! - [fn@ellippi]: Complete elliptic integral of the third kind.
//! - [fn@ellipd]: Complete elliptic integral of Legendre's type.
//! ## Legendre's incomplete integrals
//! - [fn@ellipf]: Incomplete elliptic integral of the first kind.
//! - [fn@ellipeinc]: Incomplete elliptic integral of the second kind.
//! - [fn@ellippiinc]: Incomplete elliptic integral of the third kind.
//! - [fn@ellipdinc]: Incomplete elliptic integral of Legendre's type.
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

use num_lazy::declare_nums;
declare_nums! {T}

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

// Bulirsch's integrals
pub mod bulirsch;
pub use bulirsch::cel;
pub use bulirsch::cel1;
pub use bulirsch::el1;

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

#[cfg(test)]
mod test_util;
