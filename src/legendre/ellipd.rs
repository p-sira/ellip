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

use crate::{crate_util::check_nan, elliprd, StrErr};

/// Computes [complete elliptic integral of Legendre's type](https://dlmf.nist.gov/19.2.E8).
/// ```text
///           π/2
///          ⌠       sin²θ dθ
/// D(m)  =  │  _________________
///          │     _____________
///          ⌡   \╱ 1 - m sin²θ
///         0
/// ```
///
/// ## Parameters
/// - m: elliptic parameter. m ∈ ℝ, m < 1.
///
/// The elliptic modulus (k) is also frequently used instead of the parameter (m), where k² = m.
///
/// ## Domain
/// - Returns error if m > 1.
///
/// ## Graph
/// ![Complete Elliptic Integral of Legendre's Type](https://github.com/p-sira/ellip/blob/main/figures/ellipd_plot.svg?raw=true)
///
/// [Interactive Plot](https://github.com/p-sira/ellip/blob/main/figures/ellipd_plot.html)
///
/// ## Special Cases
/// - D(0) = π/4
/// - D(1) = ∞
/// - D(-∞) = 0
///
/// # Related Functions
/// - [ellipd](crate::ellipd)(m) = ([ellipk](crate::ellipk)(m) - [ellipe](crate::ellipe)(m)) / m
/// - [ellipd](crate::ellipd)(m) = [elliprd](crate::elliprd)(0, 1 - m, 1) / 3
/// - [ellipdinc](crate::ellipdinc)(π/2, m) = [ellipd](crate::ellipd)(m)
///
/// # Examples
/// ```
/// use ellip::{ellipd, util::assert_close};
///
/// assert_close(ellipd(0.5).unwrap(), 1.0068615925073927, 1e-15);
/// ```
///
/// # References
/// - Maddock, John, Paul Bristow, Hubert Holin, and Xiaogang Zhang. “Boost Math Library: Special Functions - Elliptic Integrals.” Accessed April 17, 2025. <https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/ellint.html>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn ellipd<T: Float>(m: T) -> Result<T, StrErr> {
    check_nan!(ellipd, [m]);

    if m > 1.0 {
        return Err("ellipd: m must be less than 1.");
    }

    // D evaluates to inf at m=1.
    if m == 1.0 {
        return Ok(inf!());
    }

    if m.abs() <= epsilon!() {
        return Ok(pi!() / 4.0);
    }

    // m -> -inf, sqrt(1+m sin²θ) -> m sinθ, θ=[0,pi/2], given phi = pi/2
    // then D(m) -> 1/sqrt(-m)
    if m <= min_val!() {
        return Ok(1.0 / (-m).sqrt());
    }

    Ok(elliprd(0.0, 1.0 - m, 1.0)? / 3.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compare_test_data_boost;

    fn ellipd_k(inp: &[f64]) -> f64 {
        ellipd(inp[0] * inp[0]).unwrap()
    }

    #[test]
    fn test_ellipd() {
        compare_test_data_boost!("ellipd_data.txt", ellipd_k, 2.9e-16);
    }

    #[test]
    fn test_ellipd_special_cases() {
        use std::f64::{consts::FRAC_PI_4, INFINITY, MAX, NAN, NEG_INFINITY};
        // m = 0: D(0) = pi/4
        assert_eq!(ellipd(0.0).unwrap(), FRAC_PI_4);
        // m = 1: D(1) = inf
        assert_eq!(ellipd(1.0).unwrap(), INFINITY);
        // m < 0: should be valid
        assert!(ellipd(-1.0).unwrap().is_finite());
        // m > 1: should return Err
        assert!(ellipd(1.1).is_err());
        // m = NaN: should return Err
        assert!(ellipd(NAN).is_err());
        // m = inf: should return Err
        assert!(ellipd(INFINITY).is_err());
        // m -> -inf: D(m) = 1/sqrt(-m)
        assert_eq!(ellipd(-MAX).unwrap(), 1.0 / MAX.sqrt());
        // m = -inf: D(-inf) = 0
        assert_eq!(ellipd(NEG_INFINITY).unwrap(), 0.0);
    }
}
