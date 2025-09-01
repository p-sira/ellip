/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

//! Elliptic integral functions in Bulirsch's form.

mod cel;
mod constants;
pub(crate) mod el;

pub use cel::{cel, cel1, cel2};
pub use cel::{cel1_with_const, cel2_with_const, cel_with_const};
pub use el::{el1, el2, el3};
pub use el::{el1_with_const, el2_with_const, el3_with_const};

pub use constants::BulirschConst;

#[cfg(feature = "unstable")]
pub use constants::{DefaultPrecision, HalfPrecision};
#[cfg(feature = "unstable")]
pub use el::{el1_unchecked, el2_unchecked};
