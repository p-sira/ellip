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

pub use cel::{cel, cel1, cel2};
pub use el::{el1, el2, el3};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bulirsch_const() {
        assert_eq!(<f32 as BulirschConst>::D, 7);
        assert_eq!(<f32 as BulirschConst>::ca(), 1e-3);
        assert_eq!(<f32 as BulirschConst>::cb(), 1e-9);
        assert_eq!(<f64 as BulirschConst>::D, 16);
        assert_eq!(<f64 as BulirschConst>::ca(), 1e-8);
        assert_eq!(<f64 as BulirschConst>::cb(), 1e-18);
    }
}