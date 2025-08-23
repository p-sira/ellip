/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

//! Elliptic integral functions in Carlson's form.

mod elliprc;
mod elliprd;
mod elliprf;
mod elliprg;
mod elliprj;

pub use elliprc::elliprc;
pub use elliprd::elliprd;
pub use elliprf::elliprf;
pub use elliprg::elliprg;
pub use elliprj::elliprj;

#[cfg(not(feature = "unstable"))]
pub(crate) use {
    elliprc::elliprc_unchecked, elliprd::elliprd_unchecked, elliprf::elliprf_unchecked,
    elliprg::elliprg_unchecked, elliprj::elliprj_unchecked,
};

#[cfg(feature = "unstable")]
pub use {
    elliprc::elliprc_unchecked, elliprd::elliprd_unchecked, elliprf::elliprf_unchecked,
    elliprg::elliprg_unchecked, elliprj::elliprj_unchecked,
};
