/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2026 Sira Pornsiriprasert <code@psira.me>
 */

use num_traits::Float;

use crate::{
    crate_util::check,
    legendre::{
        ellipe::_ellipe,
        ellipk::{_ellipk, ellipk_precise_unchecked},
    },
    StrErr,
};

/// Computes [complete elliptic integrals of the first and second kind](https://dlmf.nist.gov/19.2.E8) simultaneously.
///
/// Returns `(K(m), E(m))`.
///
/// ## Parameters
/// - `m`: elliptic parameter. `m ∈ ℝ`, `m ≤ 1`.
///
/// ## Domain
/// - Returns error if `m > 1`.
/// - Returns error if `m` is NaN.
///
/// ## Special Cases
/// - `ellipke(0) = (π/2, π/2)`
/// - `ellipke(1) = (∞, 1)`
/// - `ellipke(-∞) = (0, ∞)`
///
/// # Related Functions
/// - [ellipk](crate::ellipk)
/// - [ellipe](crate::ellipe)
///
/// # Examples
/// ```
/// use ellip::{ellipke, util::assert_close};
///
/// let (k, e) = ellipke(0.5).unwrap();
/// assert_close(k, 1.8540746773013719, 1e-15);
/// assert_close(e, 1.3506438810476755, 1e-15);
/// ```
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn ellipke<T: Float>(m: T) -> Result<(T, T), StrErr> {
    check!(@nan, ellipke, [m]);

    if m > 1.0 {
        return Err("ellipke: m must not be greater than 1.");
    }
    if m == 1.0 {
        return Ok((inf!(), 1.0));
    }
    if m == neg_inf!() {
        return Ok((0.0, inf!()));
    }

    if m < 0.0 {
        if m < -1e3 {
            let k = ellipk_precise_unchecked(m);
            let c = (1.0 - m).sqrt();
            let e = c * _ellipe(m / (m - 1.0))?;
            return Ok((k, e));
        }
        let c = (1.0 - m).sqrt();
        let m1 = m / (m - 1.0);
        let (k1, e1) = _ellipke(m1)?;
        return Ok((k1 / c, c * e1));
    }

    _ellipke(m)
}

#[inline]
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
fn _ellipke<T: Float>(m: T) -> Result<(T, T), StrErr> {
    if m >= 0.9 && m < 1.0 {
        return Ok(ellipke_agm(m));
    }

    Ok((_ellipk(m)?, _ellipe(m)?))
}

/// Dedicated fused AGM kernel for near-one evaluations (`0.9 <= m < 1.0`).
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub(crate) fn ellipke_agm<T: Float>(m: T) -> (T, T) {
    let mut xn = T::one();
    let mut yn = (1.0 - m).sqrt();
    let x0 = xn;
    let y0 = yn;
    let mut sum = 0.0;
    let mut sum_pow = 0.25;

    let tol = if core::mem::size_of::<T>() <= 4 {
        1e-3
    } else {
        1e-7
    };

    for _ in 0..32 {
        let diff = (xn - yn).abs();
        if diff >= tol * xn.abs() {
            let t = (xn * yn).sqrt();
            xn = (xn + yn) / 2.0;
            yn = t;
            sum_pow = sum_pow * 2.0;
            sum = sum + sum_pow * (xn - yn) * (xn - yn);
            continue;
        }
        break;
    }

    let diff = xn - yn;
    let s = xn + yn;
    let corr = 1.0 - (diff * diff) / (4.0 * s * s);
    let k = pi!() / (s * corr);
    let e = ((x0 + y0) * (x0 + y0) / 4.0 - sum) * (pi!() / s);
    (k, e)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ellipe, ellipk, test_util::linspace, util::assert_close};

    #[test]
    fn test_ellipke_special_cases() {
        assert_eq!(ellipke(1.0).unwrap(), (f64::INFINITY, 1.0));
        assert_eq!(ellipke(f64::NEG_INFINITY).unwrap(), (0.0, f64::INFINITY));
        assert!(ellipke(f64::NAN).is_err());
        assert!(ellipke(1.1).is_err());

        let (k0, e0) = ellipke(0.0).unwrap();
        assert_close(k0, std::f64::consts::FRAC_PI_2, 1e-15);
        assert_close(e0, std::f64::consts::FRAC_PI_2, 1e-15);
    }

    #[test]
    fn test_ellipke_matches_individual() {
        // Test positive range [0, 0.999]
        for m in linspace(0.0, 0.999, 100) {
            let (k, e) = ellipke(m).unwrap();
            let k_ref = ellipk(m).unwrap();
            let e_ref = ellipe(m).unwrap();
            assert_close(k, k_ref, 1e-14);
            assert_close(e, e_ref, 1e-14);
        }

        // Test negative range [-10.0, 0.0]
        for m in linspace(-10.0, -0.001, 100) {
            let (k, e) = ellipke(m).unwrap();
            let k_ref = ellipk(m).unwrap();
            let e_ref = ellipe(m).unwrap();
            assert_close(k, k_ref, 1e-14);
            assert_close(e, e_ref, 1e-14);
        }

        // Test extreme negative range
        for &m in &[-1e4, -1e8, -1e12, -1e16, -1e50] {
            let (k, e) = ellipke(m).unwrap();
            let k_ref = ellipk(m).unwrap();
            let e_ref = ellipe(m).unwrap();
            assert_close(k, k_ref, 1e-13);
            assert_close(e, e_ref, 1e-13);
        }
    }

    #[test]
    fn test_ellipke_wolfram() {
        use crate::compare_test_data_wolfram;

        fn ellipke_k(m: f64) -> Result<f64, StrErr> {
            Ok(ellipke(m)?.0)
        }
        fn ellipke_e(m: f64) -> Result<f64, StrErr> {
            Ok(ellipke(m)?.1)
        }

        compare_test_data_wolfram!("ellipk_data.csv", ellipke_k, 1, 5e-15);
        compare_test_data_wolfram!("ellipe_data.csv", ellipke_e, 1, 5e-15);
        compare_test_data_wolfram!("ellipk_neg.csv", ellipke_k, 1, 5e-15);
        compare_test_data_wolfram!("ellipe_neg.csv", ellipke_e, 1, 5e-15);
    }
}
