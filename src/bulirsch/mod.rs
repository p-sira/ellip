/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

//! Elliptic integral functions in Bulirsch's form.

mod cel;
mod constants;
mod el;

pub use cel::{cel, cel1, cel2};
pub use el::{el1, el2, el3};

#[cfg(feature = "unstable")]
pub use cel::{_cel, _cel1, _cel2};
#[cfg(feature = "unstable")]
pub use constants::{BulirschConst, DefaultPrecision, HalfPrecision};
#[cfg(feature = "unstable")]
pub use el::{_el1, _el2, _el3, el1_unchecked, el2_unchecked};
