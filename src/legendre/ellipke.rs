/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2026 Sira Pornsiriprasert <code@psira.me>
 */

use num_traits::Float;

use super::coeffs::{self, *};
use crate::{crate_util::check, polyeval, StrErr};

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
        let c = (1.0 - m).sqrt();
        let m1 = m / (m - 1.0);
        if m1 >= 0.9 {
            // Complementary parameter mc1 = 1 - m1 = 1 / (1 - m) = 1 / c^2.
            // Passing y0 = 1 / c directly avoids catastrophic cancellation near 1.
            let (k1, e1) = ellipke_agm_from_y0(1.0 / c);
            return Ok((k1 / c, c * e1));
        }
        let (k1, e1) = _ellipke(m1)?;
        return Ok((k1 / c, c * e1));
    }

    _ellipke(m)
}

#[inline]
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
fn _ellipke<T: Float>(m: T) -> Result<(T, T), StrErr> {
    match (m * 20.0).to_i64() {
        Some(0) | Some(1) => {
            let delta = m - 0.05;
            Ok((
                polyeval(delta, &coeffs::to_t(K_0_1)),
                polyeval(delta, &coeffs::to_t(E_0_1)),
            ))
        }
        Some(2) | Some(3) => {
            let delta = m - 0.15;
            Ok((
                polyeval(delta, &coeffs::to_t(K_2_3)),
                polyeval(delta, &coeffs::to_t(E_2_3)),
            ))
        }
        Some(4) | Some(5) => {
            let delta = m - 0.25;
            Ok((
                polyeval(delta, &coeffs::to_t(K_4_5)),
                polyeval(delta, &coeffs::to_t(E_4_5)),
            ))
        }
        Some(6) | Some(7) => {
            let delta = m - 0.35;
            Ok((
                polyeval(delta, &coeffs::to_t(K_6_7)),
                polyeval(delta, &coeffs::to_t(E_6_7)),
            ))
        }
        Some(8) | Some(9) => {
            let delta = m - 0.45;
            Ok((
                polyeval(delta, &coeffs::to_t(K_8_9)),
                polyeval(delta, &coeffs::to_t(E_8_9)),
            ))
        }
        Some(10) | Some(11) => {
            let delta = m - 0.55;
            Ok((
                polyeval(delta, &coeffs::to_t(K_10_11)),
                polyeval(delta, &coeffs::to_t(E_10_11)),
            ))
        }
        Some(12) | Some(13) => {
            let delta = m - 0.65;
            Ok((
                polyeval(delta, &coeffs::to_t(K_12_13)),
                polyeval(delta, &coeffs::to_t(E_12_13)),
            ))
        }
        Some(14) | Some(15) => {
            let delta = m - 0.75;
            Ok((
                polyeval(delta, &coeffs::to_t(K_14_15)),
                polyeval(delta, &coeffs::to_t(E_14_15)),
            ))
        }
        Some(16) => {
            let delta = m - 0.825;
            Ok((
                polyeval(delta, &coeffs::to_t(K_16)),
                polyeval(delta, &coeffs::to_t(E_16)),
            ))
        }
        Some(17) => {
            let delta = m - 0.875;
            Ok((
                polyeval(delta, &coeffs::to_t(K_17)),
                polyeval(delta, &coeffs::to_t(E_17)),
            ))
        }
        Some(_) => {
            if m < 1.0 {
                Ok(ellipke_agm(m))
            } else if m == 1.0 {
                Ok((inf!(), 1.0))
            } else {
                Err("ellipke: m must not be greater than 1.")
            }
        }
        None => {
            check!(@nan, ellipke, [m]);
            if m == neg_inf!() {
                return Ok((0.0, inf!()));
            }
            #[cfg(not(feature = "test_force_fail"))]
            if m > 1.0 {
                return Err("ellipke: m must not be greater than 1.");
            }
            Err("ellipke: Unexpected error.")
        }
    }
}

/// Dedicated fused AGM kernel initialized from complementary parameter y0 = sqrt(1 - m).
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub(crate) fn ellipke_agm_from_y0<T: Float>(y0: T) -> (T, T) {
    let x0 = T::one();
    let mut xn = x0;
    let mut yn = y0;
    let mut sum = 0.0;
    let mut sum_pow = 0.25;

    let tol = T::epsilon().sqrt();

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

/// Dedicated fused AGM kernel for near-one evaluations (`0.9 <= m < 1.0`).
#[inline]
pub(crate) fn ellipke_agm<T: Float>(m: T) -> (T, T) {
    ellipke_agm_from_y0((T::one() - m).sqrt())
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

        // m < 1.0 (valid path falling through to ellipke_agm)
        let (k, e) = ellipke(0.95_f64).unwrap();
        assert!(k > 0.0 && e > 0.0);
        // m == 1.0
        assert_eq!(ellipke(1.0_f64).unwrap(), (f64::INFINITY, 1.0));
        // m > 1.0 (finite number fitting in i64 bounds)
        assert_eq!(
            ellipke(1.1_f64).unwrap_err(),
            "ellipke: m must not be greater than 1."
        );
        // m is NaN (handled by check! macro)
        assert_eq!(
            ellipke(f64::NAN).unwrap_err().to_lowercase(),
            "ellipke: arguments cannot be nan."
        );
        // m == -inf
        assert_eq!(ellipke(f64::NEG_INFINITY).unwrap(), (0.0, f64::INFINITY));
        // m > 1.0
        assert_eq!(
            ellipke(f64::INFINITY).unwrap_err(),
            "ellipke: m must not be greater than 1."
        );
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

    #[test]
    fn test_ellipke_f32_negative_accuracy() {
        let (k_f64, e_f64) = ellipke(-1000.0_f64).unwrap();
        let (k_f32, e_f32) = ellipke(-1000.0_f32).unwrap();

        let rel_err_k = ((k_f32 as f64) - k_f64).abs() / k_f64;
        let rel_err_e = ((e_f32 as f64) - e_f64).abs() / e_f64;

        // Previously at m = -1000, rel_err_k was ~2.94e-6 due to cancellation near 1.
        // With complementary parameter passed directly, relative error is within single precision (~2e-7).
        assert!(
            rel_err_k < 5e-7,
            "K relative error in f32 too high: {}",
            rel_err_k
        );
        assert!(
            rel_err_e < 5e-7,
            "E relative error in f32 too high: {}",
            rel_err_e
        );
    }

    #[test]
    fn test_ellipke_boundaries_and_extremes() {
        // Branch boundaries: m = 0.0, 0.9 (poly/AGM boundary), -9.0 (negative poly/AGM boundary), 1.0
        for &m in &[
            -1e30, -1e20, -1e10, -1000.0, -9.0, -8.999, -1.0, -0.001, 0.0, 0.5, 0.8999, 0.9, 0.95,
            0.9999,
        ] {
            let (k, e) = ellipke(m).unwrap();
            let k_ref = ellipk(m).unwrap();
            let e_ref = ellipe(m).unwrap();
            assert_close(k, k_ref, 1e-12);
            assert_close(e, e_ref, 1e-12);
        }

        // f32 extremes and boundaries
        for &m in &[
            -1e30_f32,
            -1e10_f32,
            -1000.0_f32,
            -9.0_f32,
            -1.0_f32,
            0.0_f32,
            0.5_f32,
            0.9_f32,
            0.999_f32,
        ] {
            let (k, e) = ellipke(m).unwrap();
            assert!(k.is_finite() && k > 0.0);
            assert!(e.is_finite() && e > 0.0);
        }
    }
}

#[cfg(feature = "test_force_fail")]
crate::test_force_unreachable! {
    assert_eq!(_ellipke(f64::INFINITY), Err("ellipke: Unexpected error."));
}
