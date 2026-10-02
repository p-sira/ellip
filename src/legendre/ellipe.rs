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
use crate::{crate_util::check, polyeval, StrErr};

/// Computes [complete elliptic integral of the second kind](https://dlmf.nist.gov/19.2.E8).
/// ```text
///           π/2
///          ⌠     ___________
/// E(m)  =  │ \╱ 1 - m sin²θ  dθ
///          ⌡
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
/// ![Complete Elliptic Integral of the Second Kind](https://github.com/p-sira/ellip/blob/main/figures/ellipe.svg?raw=true)
///
/// [Interactive Plot](https://p-sira.github.io/ellippy/_static/figures/ellipe.html)
///
/// ## Special Cases
/// - E(0) = π/2
/// - E(1) = 1
/// - E(-∞) = ∞
///
/// # Related Functions
/// - [ellipe](crate::ellipe)(m) = 2 [elliprg](crate::elliprg)(0, 1 - m, 1)
/// - [ellipeinc](crate::ellipeinc)(π/2, m) = [ellipe](crate::ellipe)(m)
///
/// # Examples
/// ```
/// use ellip::{ellipe, util::assert_close};
///
/// assert_close(ellipe(0.5).unwrap(), 1.3506438810476755, 1e-15);
/// ```
///
/// # References
/// - Maddock, John, Paul Bristow, Hubert Holin, and Xiaogang Zhang. “Boost Math Library: Special Functions - Elliptic Integrals.” Accessed April 17, 2025. <https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/ellint.html>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
/// - Abramowitz, Milton, and Irene A. Stegun. Handbook of Mathematical Functions: With Formulas, Graphs and Mathematical Tables. Unabridged, Unaltered and corr. Republ. of the 1964 ed. With Conference on mathematical tables, National science foundation, and Massachusetts institute of technology. Dover Books on Advanced Mathematics. Dover publ, 1972.
/// - The SciPy community. “Scipy.Special.Ellipe — SciPy v1.16.0 Manual.” Accessed July 28, 2025. <https://docs.scipy.org/doc/scipy-1.16.0/reference/generated/scipy.special.ellipe.html>.
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn ellipe<T: Float>(m: T) -> Result<T, StrErr> {
    let mut m = m;
    let mut c = 1.0;
    if m < 0.0 {
        if m == neg_inf!() {
            return Ok(inf!());
        }
        // Negative m: Abramowitz & Stegun, 1972
        c = c * (1.0 - m).sqrt();
        m = m / (m - 1.0);
    }

    Ok(c * _ellipe(m)?)
}

#[inline]
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub(crate) fn _ellipe<T: Float>(m: T) -> Result<T, StrErr> {
    match (m * 20.0).to_i64() {
        Some(0) | Some(1) => Ok(polyeval(m - 0.05, &coeffs::to_t(E_0_1))),
        Some(2) | Some(3) => Ok(polyeval(m - 0.15, &coeffs::to_t(E_2_3))),
        Some(4) | Some(5) => Ok(polyeval(m - 0.25, &coeffs::to_t(E_4_5))),
        Some(6) | Some(7) => Ok(polyeval(m - 0.35, &coeffs::to_t(E_6_7))),
        Some(8) | Some(9) => Ok(polyeval(m - 0.45, &coeffs::to_t(E_8_9))),
        Some(10) | Some(11) => Ok(polyeval(m - 0.55, &coeffs::to_t(E_10_11))),
        Some(12) | Some(13) => Ok(polyeval(m - 0.65, &coeffs::to_t(E_12_13))),
        Some(14) | Some(15) => Ok(polyeval(m - 0.75, &coeffs::to_t(E_14_15))),
        Some(16) => Ok(polyeval(m - 0.825, &coeffs::to_t(E_16))),
        Some(17) => Ok(polyeval(m - 0.875, &coeffs::to_t(E_17))),
        Some(_) => ellipe_precise(m),
        None => {
            check!(@nan, ellipe, [m]);
            #[cfg(not(feature = "test_force_fail"))]
            if m > 1.0 {
                // Infinity cases
                return Err("ellipe: m must not be greater than 1.");
            }
            Err("ellipe: Unexpected error.")
        }
    }
}

#[inline]
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
fn ellipe_precise<T: Float>(m: T) -> Result<T, StrErr> {
    // Special cases: https://dlmf.nist.gov/19.6.E1
    if m >= 1.0 {
        if m == 1.0 {
            return Ok(1.0);
        }
        return Err("ellipe: m must not be greater than 1.");
    }

    let mut xn = T::one();
    let mut yn = (1.0 - m).sqrt();
    let x0 = xn;
    let y0 = yn;
    let mut sum = 0.0;
    let mut sum_pow = 0.25;

    while (xn - yn).abs() >= 2.7 * epsilon!() * xn.abs() {
        let t = (xn * yn).sqrt();
        xn = (xn + yn) / 2.0;
        yn = t;
        sum_pow = sum_pow * 2.0;
        sum = sum + sum_pow * (xn - yn) * (xn - yn);
    }
    let rf = pi!() / (xn + yn);
    Ok(((x0 + y0) * (x0 + y0) / 4.0 - sum) * rf)
}

#[cfg(not(feature = "test_force_fail"))]
#[cfg(test)]
mod tests {
    use core::f64;

    use super::*;
    use crate::{compare_test_data_boost, compare_test_data_wolfram};

    #[test]
    fn test_ellipe() {
        compare_test_data_boost!("ellipe_data.txt", ellipe, 1, f64::EPSILON);
    }

    #[test]
    fn test_ellipe_wolfram() {
        compare_test_data_wolfram!("./tests/data/coverage", "ellipe_cov.csv", ellipe, 1, 7e-16);
    }

    #[test]
    fn test_ellipe_special_cases() {
        use std::f64::{consts::FRAC_PI_2, INFINITY, NAN, NEG_INFINITY};
        // m > 1: should return Err
        assert_eq!(ellipe(1.1), Err("ellipe: m must not be greater than 1."));
        // m = 0: E(0) = pi/2
        assert_eq!(ellipe(0.0).unwrap(), FRAC_PI_2);
        // m = 1: E(1) = 1
        assert_eq!(ellipe(1.0).unwrap(), 1.0);
        // m < 0: should be valid, compare with reference value
        assert!(ellipe(-1.0).unwrap().is_finite());
        // m = NaN: should return Err
        assert_eq!(ellipe(NAN), Err("ellipe: Arguments cannot be NAN."));
        // m = inf: should return Err
        assert_eq!(
            ellipe(INFINITY),
            Err("ellipe: m must not be greater than 1.")
        );
        // m = -inf: E(-inf) = inf
        assert_eq!(ellipe(NEG_INFINITY).unwrap(), INFINITY);
    }
}

#[cfg(feature = "test_force_fail")]
crate::test_force_unreachable! {
    assert_eq!(ellipe(f64::INFINITY), Err("ellipe: Unexpected error."));
}
