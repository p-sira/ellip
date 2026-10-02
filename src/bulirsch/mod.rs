/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

//! Elliptic integral functions in Bulirsch's form.

mod cel;
mod cel1;
mod cel2;
mod cel3;
mod constants;
mod el1;
mod el2;
pub(crate) mod el3;

pub use cel::{cel, cel_with_const};
pub use cel1::{cel1, cel1_with_const};
pub use cel2::{cel2, cel2_with_const};
pub use cel3::{cel3, cel3_with_const};
pub use el1::{el1, el1_with_const};
pub use el2::{el2, el2_with_const};
pub use el3::{el3, el3_with_const};

pub use constants::BulirschConst;

#[cfg(feature = "unstable")]
pub use constants::{DefaultPrecision, HalfPrecision};
#[cfg(feature = "unstable")]
pub use el1::el1_unchecked;
#[cfg(feature = "unstable")]
pub use el2::el2_unchecked;

/// Maximum number of iterations for [cel](crate::cel), [cel1](crate::cel1), and [cel2](crate::cel2).
#[cfg(not(feature = "test_force_fail"))]
pub(crate) const MAX_ITERATION: i16 = 10;
/// Maximum number of iterations for [cel](crate::cel), [cel1](crate::cel1), and [cel2](crate::cel2).
#[cfg(feature = "test_force_fail")]
pub(crate) const MAX_ITERATION: i16 = 1;

/// Maximum number of iterations for [el1](crate::el1), [el2](crate::el2), and [el3](crate::el3).
#[cfg(not(feature = "test_force_fail"))]
pub(crate) const N_MAX_ITERATIONS: usize = 10;
/// Maximum number of iterations for [el1](crate::el1), [el2](crate::el2), and [el3](crate::el3).
#[cfg(feature = "test_force_fail")]
pub(crate) const N_MAX_ITERATIONS: usize = 1;
