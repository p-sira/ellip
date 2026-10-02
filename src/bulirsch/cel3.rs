/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use num_traits::Float;

use crate::{
    bulirsch::{cel1::cel1_with_const, constants::BulirschConst},
    crate_util::{case, check, declare},
    StrErr,
};

/// Computes [complete elliptic integral of the third kind in Bulirsch's form](https://dlmf.nist.gov/19.2#iii).
/// ```text
///                     π/2
///                    ⌠                     dϑ
/// cel3(kc, p)     =  ⎮ ─────────────────────────────────────────────── ⋅ dϑ
///                    ⎮                         ______________________
///                    ⌡ (cos²(ϑ) + p sin²(ϑ)) ╲╱ cos²(ϑ) + kc² sin²(ϑ)
///                   0
/// ```
///
/// ## Parameters
/// - `kc`: complementary modulus. `kc ∈ ℝ`, `kc ≠ 0`.
/// - `p ∈ ℝ`, `p ≠ 0`
///
/// ## Domain
/// - Returns error if `kc = 0` or `p = 0`.
/// - Returns the Cauchy principal value for `p < 0`.
/// - Returns error if more than one arguments are infinite.
///
/// ## Graph
/// ![Complete Elliptic Integral of the Third Kind in Bulirsch's Form](https://github.com/p-sira/ellip/blob/main/figures/cel.svg?raw=true)
///
/// ## Special Cases
/// - `cel3(kc, 1) = cel1(kc)`
/// - `cel3(kc, p) = 0` for `|kc| = ∞`
/// - `cel3(kc, p) = 0` for `|p| = ∞`
///
/// # Related Functions
/// With `kc² = 1 - m` and `p = 1 - n`,
/// - [ellippi](crate::ellippi)(n, m) = [cel](crate::cel)(kc, p, 1, 1) = [cel3](crate::cel3)(kc, p)
/// - [cel3](crate::cel3)(kc, 1) = [cel1](crate::cel1)(kc) = [ellipk](crate::ellipk)(m)
///
/// # Examples
/// ```
/// use ellip::{cel3, util::assert_close};
///
/// assert_close(cel3(0.5, 0.25).unwrap(), 4.844224110273839, 1e-15);
/// ```
pub fn cel3<T: Float>(kc: T, p: T) -> Result<T, StrErr> {
    if core::mem::size_of::<T>() <= 4 {
        cel3_with_const::<T, f32>(kc, p)
    } else {
        cel3_with_const::<T, f64>(kc, p)
    }
}

/// Computes [cel3]. Control the precision using [BulirschConst].
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub fn cel3_with_const<T: Float, C: BulirschConst<T>>(kc: T, p: T) -> Result<T, StrErr> {
    check!(@nan, cel3, [kc, p]);
    check!(@zero, cel3, [kc, p]);
    check!(@multi, cel3, "infinite", is_infinite, [kc, p]);
    case!(@any [kc.abs(), p.abs()] == inf!(), T::zero());

    if p == 1.0 {
        return cel1_with_const::<T, C>(kc);
    }

    let mut kc = kc.abs();
    let mut aa: T;
    let mut bb: T;
    let mut pp: T = p;
    declare!(mut [f, g]);

    let mut e = kc;
    let mut m = 1.0;

    if pp > 0.0 {
        aa = 1.0;
        pp = pp.sqrt();
        bb = 1.0 / pp;
    } else {
        g = 1.0 - pp;
        let kc2 = kc * kc;
        if kc2.is_finite() {
            f = kc2 - pp;
            let q = 1.0 - kc2;
            pp = (f / g).sqrt();
            aa = 0.0;
            bb = -q / (g * pp);
        } else {
            let g_sqrt = g.sqrt();
            pp = kc / g_sqrt;
            aa = 0.0;
            bb = kc / g_sqrt;
        }
    }

    let mut ans = T::nan();
    for _ in 0..MAX_CEL3_ITERATION {
        f = aa;
        let inv_pp = 1.0 / pp;
        aa = bb * inv_pp + aa;
        g = if e.is_infinite() {
            (kc * inv_pp) * m
        } else {
            e * inv_pp
        };
        bb = 2.0 * (f * g + bb);
        pp = g + pp;
        g = m;
        m = kc + m;

        if (g - kc).abs() > g * C::ca() {
            kc = 2.0
                * if e.is_infinite() {
                    kc.sqrt() * m.sqrt()
                } else {
                    e.sqrt()
                };
            e = kc * m;
            continue;
        }

        ans = pi_2!() * ((aa * m + bb) / m) / (m + pp);
        break;
    }

    if ans.is_finite() {
        return Ok(ans);
    }
    Err("cel3: Failed to converge.")
}

#[cfg(not(feature = "test_force_fail"))]
const MAX_CEL3_ITERATION: i16 = 32;
#[cfg(feature = "test_force_fail")]
const MAX_CEL3_ITERATION: i16 = 1;

#[cfg(not(feature = "test_force_fail"))]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{assert_close, cel, cel1, test_util::linspace};

    #[test]
    fn test_cel3() {
        for kc in linspace(0.05, 5.0, 20) {
            // p = 1 matches cel1
            assert_close! {cel3(kc, 1.0).unwrap(), cel1(kc).unwrap(), 1e-15};

            for p in linspace(0.05, 5.0, 20) {
                let actual = cel3(kc, p).unwrap();
                let expected = cel(kc, p, 1.0, 1.0).unwrap();
                assert_close! {actual, expected, 1e-15};
            }

            for p in linspace(-5.0, -0.05, 20) {
                let actual = cel3(kc, p).unwrap();
                let expected = cel(kc, p, 1.0, 1.0).unwrap();
                assert_close! {actual, expected, 1e-14};
            }
        }
    }

    #[test]
    fn test_cel3_special_cases() {
        use std::f64::{INFINITY, NAN, NEG_INFINITY};
        assert_eq!(cel3(0.0, 1.0), Err("cel3: kc cannot be zero."));
        assert_eq!(cel3(1.0, 0.0), Err("cel3: p cannot be zero."));
        assert_eq!(cel3(INFINITY, 1.0).unwrap(), 0.0);
        assert_eq!(cel3(NEG_INFINITY, 1.0).unwrap(), 0.0);
        assert_eq!(cel3(NAN, 1.0), Err("cel3: Arguments cannot be NAN."));

        // Both infinite: returns Err
        assert!(cel3(INFINITY, INFINITY).is_err());

        // Extreme logarithmically spaced values:
        assert!(cel3(1e-300, 0.5).is_ok());
        let val_tiny = cel3(1e-300, 0.5).unwrap();
        assert!(val_tiny.is_finite() && val_tiny > 0.0);

        // cel3(1e300, 0.5) correctly returns Failed to converge, not Ok(0.0) from mutated kc
        assert_eq!(cel3(1e300, 0.5), Err("cel3: Failed to converge."));
    }

    #[test]
    fn test_cel3_independent_reference() {
        use crate::{ellipk, ellippi};

        // For kc in (0, 1), kc^2 = 1 - m => m = 1 - kc^2.
        // With n = 1 - p: cel3(kc, p) = ellippi(1 - p, 1 - kc^2).
        // ellippi uses Carlson RF and RJ, an independent recurrence.
        for kc in linspace(0.1, 0.9, 10) {
            let m = 1.0 - kc * kc;
            // p = 1 matches ellipk(m)
            assert_close! {cel3(kc, 1.0).unwrap(), ellipk(m).unwrap(), 1e-14};

            // Positive p (n < 1 and n > 1):
            for p in [0.2, 0.5, 0.8, 1.5, 2.0] {
                let actual = cel3(kc, p).unwrap();
                let expected = ellippi(1.0 - p, m).unwrap();
                assert_close! {actual, expected, 1e-12};
            }
        }
    }

    #[test]
    fn test_cel3_f32_and_extremes() {
        // f32 evaluation
        for &kc in &[0.1_f32, 0.5_f32, 1.0_f32, 2.0_f32] {
            for &p in &[-2.0_f32, -0.5_f32, 0.25_f32, 1.0_f32, 3.0_f32] {
                let val_f32 = cel3(kc, p).unwrap();
                let val_f64 = cel3(kc as f64, p as f64).unwrap();
                assert_close! {val_f32 as f64, val_f64, 1e-4};
            }
        }

        // Logarithmically spaced extremes
        for &kc in &[1e-300, 1e-100, 1e-50, 1e-10, 1e-1, 1.0, 10.0, 100.0] {
            let res = cel3(kc, 0.5);
            assert!(res.is_ok());
            let val = res.unwrap();
            assert!(val.is_finite() && val > 0.0);
        }

        // Negative p extremes
        for &p in &[-1e5, -100.0, -10.0, -0.5, -0.01] {
            let res = cel3(0.5, p);
            assert!(res.is_ok());
            assert!(res.unwrap().is_finite());
        }
    }
}

#[cfg(feature = "test_force_fail")]
crate::test_force_unreachable! {
    assert_eq!(cel3(1e300, 0.5), Err("cel3: Failed to converge."));
}
