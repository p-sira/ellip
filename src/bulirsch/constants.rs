/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */
use num_traits::Float;

pub trait BulirschConst<T: Float> {
    /// Number of significant figures.
    ///
    /// **D** is `7` for `f32` and `16` for `f64`.
    #[allow(dead_code)]
    const D: i32;

    /// D-2
    const ND: usize;

    /// 1e(-D/2)
    fn ca() -> T;

    /// 1e(-(D+2))
    fn cb() -> T;

    /// Limit of the kc and p to use ellippiinc
    /// defaults to 1e-4 * CA
    fn lim_kc_p() -> T;
}

macro_rules! impl_bulirsch_const {
    ($d:literal, $ca:literal, $cb:literal, $lim_kc_p:literal) => {
        const D: i32 = $d;
        const ND: usize = $d - 2;
        fn ca() -> T {
            T::from($ca).unwrap()
        }
        fn cb() -> T {
            T::from($cb).unwrap()
        }
        fn lim_kc_p() -> T {
            T::from($lim_kc_p).unwrap()
        }
    };
    (@type $type:ty, {D: $d:literal, CA: $ca:literal, CB: $cb:literal, LIM: $lim_kc_p:literal}) => {
        impl<T: Float> BulirschConst<T> for $type {
            impl_bulirsch_const!($d, $ca, $cb, $lim_kc_p);
        }
    };
    ($struct:ident, {D: $d:literal, CA: $ca:literal, CB: $cb:literal, LIM: $lim_kc_p:literal}) => {
        pub struct $struct;
        impl<T: Float> BulirschConst<T> for $struct {
            impl_bulirsch_const!($d, $ca, $cb, $lim_kc_p);
        }
    };
}

impl_bulirsch_const!(@type f32, {D: 7, CA: 1e-3, CB: 1e-9, LIM: 1e-7});
impl_bulirsch_const!(@type f64, {D: 16, CA: 1e-8, CB: 1e-18, LIM: 1e-12});
#[cfg(feature = "unstable")]
impl_bulirsch_const!(HalfPrecision, {D: 7, CA: 1e-3, CB: 1e-9, LIM: 1e-7});
impl_bulirsch_const!(DefaultPrecision, {D: 16, CA: 1e-8, CB: 1e-18, LIM: 1e-12});

#[cfg(not(feature = "test_force_fail"))]
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bulirsch_const() {
        assert_eq!(<f32 as BulirschConst<f32>>::D, 7);
        assert_eq!(<f32 as BulirschConst<f32>>::ca(), 1e-3);
        assert_eq!(<f32 as BulirschConst<f32>>::cb(), 1e-9);
        assert_eq!(<f32 as BulirschConst<f32>>::lim_kc_p(), 1e-7);
        assert_eq!(<f64 as BulirschConst<f64>>::D, 16);
        assert_eq!(<f64 as BulirschConst<f64>>::ca(), 1e-8);
        assert_eq!(<f64 as BulirschConst<f64>>::cb(), 1e-18);
        assert_eq!(<f64 as BulirschConst<f64>>::lim_kc_p(), 1e-12);
    }
}
