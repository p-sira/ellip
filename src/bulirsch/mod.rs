/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

//! Elliptic integral functions in Bulirsch's form.

pub trait BulirschConst {
    /// Number of significant figures.
    /// 
    /// **D** is `7` for `f32` and `16` for `f64`.
    const D: i32;

    /// 10^(-D/2)
    fn ca() -> Self;

    /// 10^(-(D+2))
    fn cb() -> Self;
}

impl BulirschConst for f32 {
    const D: i32 = 7;

    fn ca() -> Self {
        10_f32.powi(-Self::D / 2)
    }

    fn cb() -> Self {
        10_f32.powi(-(Self::D + 2))
    }
}

impl BulirschConst for f64 {
    const D: i32 = 16;

    fn ca() -> Self {
        10_f64.powi(-Self::D / 2)
    }

    fn cb() -> Self {
        10_f64.powi(-(Self::D + 2))
    }
}

mod cel;
mod el;

pub use cel::cel;
pub use el::el1;
