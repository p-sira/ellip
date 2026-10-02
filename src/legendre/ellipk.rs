/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Boost Math, see LICENSE in this directory.
 */

// Original header from Boost Math
//  Copyright (c) 2006 Xiaogang Zhang
//  Copyright (c) 2006 John Maddock
//  Copyright (c) 2024 Matt Borland
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0.

use num_traits::Float;

use super::coeffs::{self, *};
use crate::{
    crate_util::{check, declare},
    polyeval, StrErr,
};

/// Computes [complete elliptic integral of the first kind](https://dlmf.nist.gov/19.2.E8).
/// ```text
///           π/2
///          ⌠          dθ
/// K(m)  =  │  _________________
///          │     _____________
///          ⌡   \╱ 1 - m sin²θ
///         0
/// ```
///
/// ## Parameters
/// - m: elliptic parameter. m ∈ ℝ, m ≤ 1.
///
/// The elliptic modulus (k) is also frequently used instead of the parameter (m), where k² = m.
///
/// ## Domain
/// - Returns error if m > 1.
///
/// ## Graph
/// ![Complete Elliptic Integral of the First Kind](https://github.com/p-sira/ellip/blob/main/figures/ellipk.svg?raw=true)
///
/// [Interactive Plot](https://p-sira.github.io/ellippy/_static/figures/ellipk.html)
///
/// ## Special Cases
/// - K(0) = π/2
/// - K(1) = ∞
/// - K(-∞) = 0
///
/// # Related Functions
/// - [ellipk](crate::ellipk)(m) = [elliprf](crate::elliprf)(0, 1 - m, 1)
/// - [ellipf](crate::ellipf)(π/2, m) = [ellipk](crate::ellipk)(m)
///
/// # Examples
/// ```
/// use ellip::{ellipk, util::assert_close};
///
/// assert_close(ellipk(0.5).unwrap(), 1.8540746773013719, 1e-15);
/// ```
///
/// # References
/// - Maddock, John, Paul Bristow, Hubert Holin, and Xiaogang Zhang. “Boost Math Library: Special Functions - Elliptic Integrals.” Accessed April 17, 2025. <https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/ellint.html>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn ellipk<T: Float>(m: T) -> Result<T, StrErr> {
    // Negative parameters use the AGM, never the polynomial selector.
    if m < 0.0 && m.is_finite() {
        return ellipk_precise(m);
    }
    _ellipk(m)
}

#[inline]
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub(crate) fn _ellipk<T: Float>(m: T) -> Result<T, StrErr> {
    match (m * 20.0).to_i64() {
        Some(0) | Some(1) => Ok(polyeval(m - 0.05, &coeffs::to_t(K_0_1))),
        Some(2) | Some(3) => Ok(polyeval(m - 0.15, &coeffs::to_t(K_2_3))),
        Some(4) | Some(5) => Ok(polyeval(m - 0.25, &coeffs::to_t(K_4_5))),
        Some(6) | Some(7) => Ok(polyeval(m - 0.35, &coeffs::to_t(K_6_7))),
        Some(8) | Some(9) => Ok(polyeval(m - 0.45, &coeffs::to_t(K_8_9))),
        Some(10) | Some(11) => Ok(polyeval(m - 0.55, &coeffs::to_t(K_10_11))),
        Some(12) | Some(13) => Ok(polyeval(m - 0.65, &coeffs::to_t(K_12_13))),
        Some(14) | Some(15) => Ok(polyeval(m - 0.75, &coeffs::to_t(K_14_15))),
        Some(16) => Ok(polyeval(m - 0.825, &coeffs::to_t(K_16))),
        Some(17) => Ok(polyeval(m - 0.875, &coeffs::to_t(K_17))),
        Some(_) => ellipk_precise(m),
        None => {
            check!(@nan, ellipk, [m]);
            if m == neg_inf!() {
                return Ok(0.0);
            }
            #[cfg(not(feature = "test_force_fail"))]
            if m > 1.0 {
                // Also handles inf
                return Err("ellipk: m must not be greater than 1.");
            }
            Err("ellipk: Unexpected error.")
        }
    }
}

#[inline]
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub(crate) fn ellipk_precise<T: Float>(m: T) -> Result<T, StrErr> {
    // Special cases: https://dlmf.nist.gov/19.6.E1
    if m >= 1.0 {
        if m == 1.0 {
            return Ok(inf!());
        }
        return Err("ellipk: m must not be greater than 1.");
    }

    Ok(ellipk_precise_unchecked(m))
}

/// Based on elliprf(1, 1-m, 0)
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub fn ellipk_precise_unchecked<T: Float>(m: T) -> T {
    declare!(mut [xn = T::one(), yn = (T::one() - m).sqrt(), t]);

    let tol = if core::mem::size_of::<T>() <= 4 {
        1e-3
    } else {
        1e-7
    };

    for _ in 0..MAX_ITERATION {
        let diff = (xn - yn).abs();
        if diff >= tol * xn.abs() {
            t = (xn * yn).sqrt();
            xn = (xn + yn) / 2.0;
            yn = t;
            continue;
        }
        break;
    }

    let diff = xn - yn;
    let sum = xn + yn;
    let corr = 1.0 - (diff * diff) / (4.0 * sum * sum);
    pi!() / (sum * corr)
}

const MAX_ITERATION: usize = 32;

#[cfg(not(feature = "test_force_fail"))]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::compare_test_data_boost;
    use crate::compare_test_data_wolfram;

    #[test]
    fn test_ellipk_boost() {
        compare_test_data_boost!("ellipk_data.txt", ellipk, 1, 2.8e-16);
    }

    #[test]
    fn test_ellipk_wolfram() {
        compare_test_data_wolfram!("./tests/data/coverage", "ellipk_cov.csv", ellipk, 1, 2e-12);
        compare_test_data_wolfram!("ellipk_data.csv", ellipk, 1, 5e-15);
    }

    #[test]
    fn test_ellipk_special_cases() {
        use std::f64::{consts::FRAC_PI_2, INFINITY, NAN, NEG_INFINITY};
        // m = 0: K(0) = pi/2
        assert_eq!(ellipk(0.0).unwrap(), FRAC_PI_2);
        // m = 1: K(1) = inf
        assert_eq!(ellipk(1.0).unwrap(), INFINITY);
        // m < 0: should be valid, compare with reference value
        assert!(ellipk(-1.0).unwrap().is_finite());
        // m > 1: should return Err
        assert_eq!(ellipk(1.1), Err("ellipk: m must not be greater than 1."));
        // m = NaN: should return Err
        assert_eq!(ellipk(NAN), Err("ellipk: Arguments cannot be NAN."));
        // m = inf: should return Err
        assert_eq!(
            ellipk(INFINITY),
            Err("ellipk: m must not be greater than 1.")
        );
        // m = -inf: K(-inf) = 0
        assert_eq!(ellipk(NEG_INFINITY).unwrap(), 0.0);
    }

    // Regression for audit finding A10: https://github.com/p-sira/ellip/pull/124
    #[test]
    fn test_large_negative_parameters() {
        crate::assert_close!(ellipk(-1e18).unwrap(), 2.2109560198066302e-8, 2e-15);
        crate::assert_close!(ellipk(-1e100).unwrap(), 1.1651554901082218e-48, 3e-15);
        crate::assert_close!(ellipk(-f64::MAX).unwrap(), 2.6572401146362276e-152, 3e-15);
        assert_eq!(ellipk(f64::NEG_INFINITY).unwrap(), 0.0);
        assert!(ellipk(f64::NAN).is_err());
    }
}

#[cfg(feature = "test_force_fail")]
crate::test_force_unreachable! {
    assert_eq!(ellipk(f64::INFINITY), Err("ellipk: Unexpected error."));
}
