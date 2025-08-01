/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */
use num_traits::Float;

#[deprecated = "Change signature in the next minor version and will be moved to unstable feature."]
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

macro_rules! _impl_bulirsch_const {
    ($type: ty, $d: literal) => {
        impl BulirschConst for $type {
            const D: i32 = $d;

            fn ca() -> Self {
                (10 as $type).powi(-<Self as BulirschConst>::D / 2)
            }

            fn cb() -> Self {
                (10 as $type).powi(-(<Self as BulirschConst>::D + 2))
            }
        }
    };
}

_impl_bulirsch_const!(f32, 7);
_impl_bulirsch_const!(f64, 16);

/// Hides as internal function until next minor release
#[deprecated = "New feature. Release on next minor version as unstable feature."]
pub trait _BulirschConst<T: Float> {
    /// Number of significant figures.
    ///
    /// **D** is `7` for `f32` and `16` for `f64`.
    const D: i32;

    /// D-2
    const ND: usize;

    /// 10^(-D/2)
    fn ca() -> T;

    /// 10^(-(D+2))
    fn cb() -> T;
}

macro_rules! impl_bulirsch_const {
    ($d:literal, $ca:literal, $cb:literal) => {
        const D: i32 = $d;
        const ND: usize = $d - 2;
        fn ca() -> T {
            T::from($ca).unwrap()
        }
        fn cb() -> T {
            T::from($cb).unwrap()
        }
    };
    (@type $type:ty, {D: $d:literal, CA: $ca:literal, CB: $cb:literal}) => {
        impl<T: Float> _BulirschConst<T> for $type {
            impl_bulirsch_const!($d, $ca, $cb);
        }
    };
    ($struct:ident, {D: $d:literal, CA: $ca:literal, CB: $cb:literal}) => {
        pub struct $struct;
        impl<T: Float> _BulirschConst<T> for $struct {
            impl_bulirsch_const!($d, $ca, $cb);
        }
    };
}

impl_bulirsch_const!(@type f32, {D: 7, CA: 1e-3, CB: 1e-9});
impl_bulirsch_const!(@type f64, {D: 16, CA: 1e-8, CB: 1e-18});
impl_bulirsch_const!(HalfPrecision, {D: 7, CA: 1e-3, CB: 1e-9});
impl_bulirsch_const!(DefaultPrecision, {D: 16, CA: 1e-8, CB: 1e-18});

#[cfg(not(feature = "reduce-iteration"))]
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bulirsch_const() {
        assert_eq!(<f32 as _BulirschConst<f32>>::D, 7);
        assert_eq!(<f32 as _BulirschConst<f32>>::ca(), 1e-3);
        assert_eq!(<f32 as _BulirschConst<f32>>::cb(), 1e-9);
        assert_eq!(<f64 as _BulirschConst<f64>>::D, 16);
        assert_eq!(<f64 as _BulirschConst<f64>>::ca(), 1e-8);
        assert_eq!(<f64 as _BulirschConst<f64>>::cb(), 1e-18);
    }
}
