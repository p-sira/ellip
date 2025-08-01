/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

//! Elliptic integral functions in Bulirsch's form.

mod cel;
mod constants;
mod el;

pub use cel::{cel, cel1, cel2};
pub use constants::BulirschConst;
pub(crate) use constants::{DefaultPrecision, _BulirschConst};
pub use el::{el1, el2, el3};
