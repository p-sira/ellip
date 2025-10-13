/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

//! Functions related to elliptic integrals

#[cfg(feature = "unstable")]
pub mod jacobi_zeta;
#[cfg(not(feature = "unstable"))]
mod jacobi_zeta;
pub use jacobi_zeta::jacobi_zeta;
#[cfg(feature = "unstable")]
pub mod heuman_lambda;
#[cfg(not(feature = "unstable"))]
mod heuman_lambda;
pub use heuman_lambda::heuman_lambda;
