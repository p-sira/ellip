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
//! - [fn@cel]: General complete elliptic integral in Bulirsch's form.
//! - [fn@cel1]: Complete elliptic integral of the first kind in Bulirsch's form.
//! - [fn@el1]: Incomplete elliptic integral of the first kind in Bulirsch's form.
//! - [fn@el2]: Incomplete elliptic integral of the second kind in Bulirsch's form.
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
//! References for original implementations are:
//! - NIST Digital Library, [Chapter 19: Elliptic Integrals](https://dlmf.nist.gov/19) (Carlson, 2024).
//! - Numerical calculation of elliptic integrals and elliptic functions [I](https://link.springer.com/article/10.1007/BF01397975) (Bulirsch, 1965), [II](https://doi.org/10.1007/BF01436529) (Bulirsch, 1965), and [III](https://doi.org/10.1007/BF02165405) (Bulirsch, 1969).
//! - [Cylindrical magnets and ideal solenoids](https://doi.org/10.1119/1.3256157) (Derby and Olbert, 2010).
//!
//! Unicode-style mathematical notation are created using [Diagon](https://github.com/ArthurSonzogni/Diagon).

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
pub use bulirsch::el2;

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
