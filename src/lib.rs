/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

//! **Ellip** is a standalone elliptic integral functions for Rust.
//!
//! Ellip is based on Scipy C++ implementations, which were translated from Cephes Math Library.
//!
//! # Features
//! ## Incomplete elliptic integrals
//! - [ellipf]: Incomplete elliptic integral of the first kind.
//! - [ellipeinc]: Incomplete elliptic integral of the second kind.
//! ## Complete elliptic integrals
//! - [ellipk]: Complete elliptic integral of the first kind.
//! - [ellipe]: Complete elliptic integral of the second kind.
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

#[cfg(test)]
mod test_util;
