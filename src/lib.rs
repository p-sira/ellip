/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

//! # Ellip
//! **Ellip** is a standalone elliptic integral functions for Rust.
//!
//! # Features
//! ## Incomplete elliptic integrals
//! - [fn@ellipf]: Incomplete elliptic integral of the first kind.
//! - [fn@ellipeinc]: Incomplete elliptic integral of the second kind.
//! ## Complete elliptic integrals
//! - [fn@ellipk]: Complete elliptic integral of the first kind.
//! - [fn@ellipe]: Complete elliptic integral of the second kind.
//! ## Symmetric elliptic integrals
//! - [fn@elliprf]: Symmetric elliptic integral of the first kind.
//! - [fn@elliprd]: Degenerate elliptic integral of the third kind.
//! - [fn@elliprc]: Degenerate elliptic integral of RF.
//!
//! # Acknowledgment
//! Ellip is derived from multiple mathematic libraries. We thank
//! the opensource contributors for making mathematic libraries free for all.
//! Following are the main original works used in the development of Ellip.
//! Detailed credits are available in the source code.
//! - [Scipy](https://github.com/scipy/scipy/)
//! - [Cephes Math Library](https://netlib.org/cephes/)
//! - [Russell Lab](https://github.com/cpmech/russell)
//!

mod polyeval;
use polyeval::*;

mod ellipf;
pub use ellipf::ellipf;
mod ellipeinc;
pub use ellipeinc::ellipeinc;

mod ellipe;
pub use ellipe::ellipe;
mod ellipk;
pub use ellipk::ellipk;

mod elliprf;
pub use elliprf::elliprf;
mod elliprd;
pub use elliprd::elliprd;
mod elliprc;
pub use elliprc::elliprc;

mod constants;

#[cfg(test)]
mod test_util;
