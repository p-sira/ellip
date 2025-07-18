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
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  History:
//  XZ wrote the original of this file as part of the Google
//  Summer of Code 2006.  JM modified it to fit into the
//  Boost.Math conceptual framework better, and to ensure
//  that the code continues to work no matter how many digits
//  type T has.

use crate::{ellipd, elliprd, StrErr};
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
    let sign = if phi < 0.0 { -1.0 } else { 1.0 };

    let phi = phi.abs();

    if phi >= max_val!() {
        // Need to handle infinity as a special case:
        return Ok(inf!());
    }

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

#[cfg(test)]
mod tests {
    use core::f64;

    use super::*;
    use crate::compare_test_data_boost;

    fn ellipdinc_k(inp: &[f64]) -> f64 {
        ellipdinc(inp[0], inp[1] * inp[1]).unwrap()
    }

    #[test]
    fn test_ellipd() {
        compare_test_data_boost!("ellipdinc_data.txt", ellipdinc_k, 6.4e-16);
    }

    #[test]
    fn test_ellipdinc_err() {
        use std::f64::consts::FRAC_PI_2;
        // m sin²φ > 1
        assert!(ellipdinc(FRAC_PI_2, 1.1).is_err());
    }
}
