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

use crate::{crate_util::check, ellipd, elliprd, StrErr};
use num_traits::Float;

/// Computes [incomplete elliptic integral of Legendre's type](https://dlmf.nist.gov/19.2.E6).
/// ```text
///              φ
///             ⌠       sin²θ dθ
/// D(φ, m)  =  │  _________________
///             │     _____________
///             ⌡   \╱ 1 - m sin²θ
///            0
/// ```
///
/// ## Parameters
/// - phi: amplitude angle (φ). φ ∈ ℝ.
/// - m: elliptic parameter. m ∈ ℝ.
///
/// The elliptic modulus (k) is also frequently used instead of the parameter (m), where k² = m.
///
/// ## Domain
/// - Returns error if m sin²φ > 1.
///
/// ## Graph
/// ![Incomplete Elliptic Integral of Legendre's Type](https://github.com/p-sira/ellip/blob/main/figures/ellipdinc_plot.svg?raw=true)
///
/// [Interactive Plot](https://github.com/p-sira/ellip/blob/main/figures/ellipdinc_plot.html)
///
/// ![Incomplete Elliptic Integral of Legendre's Type (3D Plot)](https://github.com/p-sira/ellip/blob/main/figures/ellipdinc_plot_3d.png?raw=true)
///
/// [Interactive 3D Plot](https://github.com/p-sira/ellip/blob/main/figures/ellipdinc_plot_3d.html)
///
/// ## Special Cases
/// - D(0, m) = 0
/// - D(φ, 0) = sin φ
/// - D(π/2, m) = D(m)
/// - D(φ, -∞) = 0
/// - D(∞, m) = ∞
/// - D(-∞, m) = -∞
///
/// # Related Functions
/// With c = csc²φ,
/// - [ellipdinc](crate::ellipdinc)(φ, m) = ([ellipf](crate::ellipf)(φ, m) - [ellipeinc](crate::ellipeinc)(φ, m)) / m
/// - [ellipdinc](crate::ellipdinc)(φ, m) = [elliprd](crate::elliprd)(c - 1, c - m, c) / 3
///
/// # Examples
/// ```
/// use ellip::{ellipdinc, util::assert_close};
/// use std::f64::consts::FRAC_PI_4;
///
/// assert_close(ellipdinc(FRAC_PI_4, 0.5).unwrap(), 0.15566274414316758, 1e-15);
/// ```
///
/// # References
/// - Maddock, John, Paul Bristow, Hubert Holin, and Xiaogang Zhang. “Boost Math Library: Special Functions - Elliptic Integrals.” Accessed April 17, 2025. <https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/ellint.html>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn ellipdinc<T: Float>(phi: T, m: T) -> Result<T, StrErr> {
    check!(@nan, ellipdinc, [phi, m]);

    let sign = if phi < 0.0 { -1.0 } else { 1.0 };

    if phi.abs() >= max_val!() {
        // Need to handle infinity as a special case:
        return Ok(phi.signum() * inf!());
    }

    // m=-inf should be undefined but as m -> -inf
    // D -> -inf
    if m == neg_inf!() {
        return Ok(neg_inf!());
    }

    let phi = phi.abs();

    if phi > 1.0 / epsilon!() {
        // Phi is so large that phi%pi is necessarily zero (or garbage),
        // just return the second part of the duplication formula:
        return Ok(sign * 2.0 * phi * ellipd(m)? / pi!());
    }

    // Carlson's algorithm works only for |phi| <= pi/2,
    // use the integrand's periodicity to normalize phi
    //
    let mut rphi = phi % pi_2!();
    let mut mm = ((phi - rphi) / pi_2!()).round();
    let mut s = 1.0;
    if mm % 2.0 > 0.5 {
        mm = mm + 1.0;
        s = -1.0;
        rphi = pi_2!() - rphi;
    }

    let sinp = rphi.sin();
    let cosp = rphi.cos();
    let sinp2 = sinp * sinp;
    let cosp2 = cosp * cosp;

    let c = 1.0 / sinp2;
    let cm1 = cosp2 / sinp2; // c - 1

    if m * sinp2 > 1.0 {
        return Err("ellipdinc: m sin²φ must be smaller than one.");
    }

    let mut result = 0.0;
    if rphi != 0.0 {
        result = s * elliprd(cm1, c - m, c)? / 3.0;
    }
    if mm != 0.0 {
        result = result + mm * ellipd(m)?;
    }

    Ok(sign * result)
}

#[cfg(not(feature = "reduce-iteration"))]
#[cfg(test)]
mod tests {
    use core::f64;

    use super::*;
    use crate::compare_test_data_boost;

    fn ellipdinc_k(inp: &[f64]) -> f64 {
        ellipdinc(inp[0], inp[1] * inp[1]).unwrap()
    }

    #[test]
    fn test_ellipdinc() {
        compare_test_data_boost!("ellipdinc_data.txt", ellipdinc_k, 6.4e-16);
    }

    #[test]
    fn test_ellipdinc_special_cases() {
        use std::f64::{
            consts::{FRAC_PI_2, FRAC_PI_4, PI},
            INFINITY, NAN, NEG_INFINITY,
        };
        // phi = pi/2, m = 1: D(pi/2, 1) = inf
        assert_eq!(ellipdinc(FRAC_PI_2, 1.0).unwrap(), INFINITY);
        // m * sin^2(phi) >= 1: should return Err
        assert!(ellipdinc(FRAC_PI_2, 2.0).is_err());
        // phi = 0: D(0, m) = 0
        assert_eq!(ellipdinc(0.0, 0.5).unwrap(), 0.0);
        // phi = pi/2, m = 0: D(pi/2, 0) = pi/4
        assert_eq!(ellipdinc(FRAC_PI_2, 0.0).unwrap(), FRAC_PI_4);
        // m < 0: should be valid
        assert!(ellipdinc(FRAC_PI_2, -1.0).unwrap().is_finite());
        // phi > 1/epsilon: D(phi, m) = 2 * phi * D(m) / pi
        assert_eq!(
            ellipdinc(1e16, 0.4).unwrap(),
            2.0 * 1e16 * ellipd(0.4).unwrap() / PI
        );
        // phi = nan or m = nan: should return Err
        assert!(ellipdinc(NAN, 0.5).is_err());
        assert!(ellipdinc(0.5, NAN).is_err());
        // phi = inf: D(inf, m) = inf
        assert_eq!(ellipdinc(INFINITY, 0.5).unwrap(), INFINITY);
        // phi = -inf: D(-inf, m) = -inf
        assert_eq!(ellipdinc(NEG_INFINITY, 0.5).unwrap(), NEG_INFINITY);
        // m = inf: should return Err
        assert!(ellipdinc(0.5, INFINITY).is_err());
        // m = -inf: D(phi, -inf) = -inf
        assert_eq!(ellipdinc(0.5, NEG_INFINITY).unwrap(), NEG_INFINITY);
    }
}
