/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

//! Elliptic integral functions in Legendre's form.

mod ellipd;
mod ellipdinc;
mod ellipe;
mod ellipeinc;
mod ellipf;
mod ellipk;
mod ellippi;
mod ellippiinc;

pub use ellipd::ellipd;
pub use ellipdinc::ellipdinc;
pub use ellipe::ellipe;
pub use ellipeinc::ellipeinc;
pub use ellipf::ellipf;
pub use ellipk::ellipk;
pub use ellippi::ellippi;
pub use ellippiinc::{ellippiinc, ellippiinc_bulirsch};

#[cfg(feature = "unstable")]
pub use {ellipeinc::ellipeinc_unchecked, ellippi::ellippi_unchecked};
