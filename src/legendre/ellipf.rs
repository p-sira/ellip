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

use crate::{
    carlson::elliprf_unchecked,
    crate_util::{case, check},
    legendre::ellipk::ellipk_precise_unchecked,
    StrErr,
};

/// Computes [incomplete elliptic integral of the first kind](https://dlmf.nist.gov/19.2.E4).
/// ```text
///              φ
///             ⌠          dθ
/// F(φ, m)  =  │  _________________
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
/// ![Incomplete Elliptic Integral of the First Kind](https://github.com/p-sira/ellip/blob/main/figures/ellipf.svg?raw=true)
///
/// [Interactive Plot](https://p-sira.github.io/ellippy/_static/figures/ellipf.html)
///
/// ![Incomplete Elliptic Integral of the First Kind (3D Plot)](https://github.com/p-sira/ellip/blob/main/figures/ellipf_3d.svg?raw=true)
///
/// [Interactive 3D Plot](https://p-sira.github.io/ellippy/_static/figures/ellipf_3d.html)
///
/// ## Special Cases
/// - F(0, m) = 0
/// - F(φ, 0) = φ
/// - F(π/2, m) = K(m)
/// - F(φ, -∞) = 0
/// - F(φ, m) = sign(φ) ∞ for |φ| = ∞
///
/// # Related Functions
/// With c = csc²φ,
/// - [ellipf](crate::ellipf)(φ,m) = [elliprf](crate::elliprf)(c - 1, c - m, c)
///
/// # Examples
/// ```
/// use ellip::{ellipf, util::assert_close};
/// use std::f64::consts::FRAC_PI_4;
///
/// assert_close(ellipf(FRAC_PI_4, 0.5).unwrap(), 0.826017876249245, 1e-15);
/// ```
///
/// # References
/// - Maddock, John, Paul Bristow, Hubert Holin, and Xiaogang Zhang. “Boost Math Library: Special Functions - Elliptic Integrals.” Accessed April 17, 2025. <https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/ellint.html>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
/// - The MathWorks, Inc. “ellipticF.” Accessed April 21, 2025. <https://www.mathworks.com/help/symbolic/sym.ellipticf.html>.
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn ellipf<T: Float>(phi: T, m: T) -> Result<T, StrErr> {
    let sign = phi.signum();
    let phi = phi.abs();

    // Large phis
    if phi > 1.0 / epsilon!() {
        if phi >= max_val!() {
            return Ok(sign * inf!());
        }
        // Phi is so large that phi%pi is necessarily zero (or garbage),
        // just return the second part of the duplication formula:
        return Ok(sign * 2.0 * phi * ellipk_precise_unchecked(m) / pi!());
    }

    // Carlson's algorithm works only for |phi| <= pi/2,
    // use the integrand's periodicity to normalize phi
    let mut rphi = phi % pi_2!();
    let mut mm = ((phi - rphi) / pi_2!()).round();

    let mut s = 1.0;
    if mm % 2.0 > 0.5 {
        mm = mm + 1.0;
        s = -1.0;
        rphi = pi_2!() - rphi;
    }

    let sphi = rphi.sin();
    let s2p = sphi * sphi;
    if m * s2p >= 1.0 {
        return Err("ellipf: m sin²φ must be smaller than one.");
    }
    let cphi = rphi.cos();
    let c2p = cphi * cphi;
    let mut ans;

    debug_assert!(s2p > min_val!() || !s2p.is_normal());
    // Use http://dlmf.nist.gov/19.25#E5, note that
    // c-1 simplifies to cot^2(rphi) which avoids cancellation.
    // Likewise c - k^2 is the same as (c - 1) + (1 - k^2).
    //
    let c = s2p.recip();
    let c_minus_one = c2p / s2p;
    let arg2 = if m != 0.0 {
        let cross = (c / m).abs();
        if cross > 0.9 && cross < 1.1 {
            c_minus_one + 1.0 - m
        } else {
            c - m
        }
    } else {
        c
    };
    ans = s * elliprf_unchecked(c_minus_one, arg2, c);

    // It should not be possible for s2p to be less than MIN value.
    // if s2p > min_val!() {
    //     // Use http://dlmf.nist.gov/19.25#E5, note that
    //     // c-1 simplifies to cot^2(rphi) which avoids cancellation.
    //     // Likewise c - k^2 is the same as (c - 1) + (1 - k^2).
    //     //
    //     let c = 1.0 / s2p;
    //     let c_minus_one = c2p / s2p;
    //     let arg2 = if m != 0.0 {
    //         let cross = (c / m).abs();
    //         if cross > 0.9 && cross < 1.1 {
    //             c_minus_one + 1.0 - m
    //         } else {
    //             c - m
    //         }
    //     } else {
    //         c
    //     };
    //     result = s * elliprf(c_minus_one, arg2, c)?;
    // } else {
    //     result = s * rphi.sin();
    // }

    if mm != 0.0 {
        ans = ans + mm * ellipk_precise_unchecked(m);
    }

    ans = sign * ans;
    if ans.is_finite() {
        #[cfg(not(feature = "test_force_fail"))]
        return Ok(ans);
    }
    check!(@nan, ellipf, [phi, m]);
    case!(m == neg_inf!(), T::zero());
    Err("ellipf: Unexpected error.")
}

#[cfg(not(feature = "test_force_fail"))]
#[cfg(all(test, not(feature = "no_std")))]
mod tests {
    use super::*;
    use crate::compare_test_data_boost;

    #[test]
    fn test_ellipf() {
        compare_test_data_boost!("ellipf_data.txt", ellipf, 2, 5.1e-16);
    }

    #[test]
    fn test_ellipf_special_cases() {
        use crate::ellipk;
        use std::f64::{
            consts::{FRAC_PI_2, PI},
            INFINITY, NAN, NEG_INFINITY,
        };
        // m * sin^2(phi) >= 1: should return Err
        assert_eq!(
            ellipf(FRAC_PI_2, 1.0),
            Err("ellipf: m sin²φ must be smaller than one.")
        );
        assert_eq!(
            ellipf(FRAC_PI_2, 2.0),
            Err("ellipf: m sin²φ must be smaller than one.")
        );
        // phi = 0: F(0, m) = 0
        assert_eq!(ellipf(0.0, 0.5).unwrap(), 0.0);
        // phi = pi/2, m = 0: F(pi/2, 0) = pi/2
        assert_eq!(ellipf(FRAC_PI_2, 0.0).unwrap(), FRAC_PI_2);
        // m < 0: should be valid
        assert!(ellipf(FRAC_PI_2, -1.0).unwrap().is_finite());
        // phi = nan or m = nan: should return Err
        assert_eq!(ellipf(NAN, 0.5), Err("ellipf: Arguments cannot be NAN."));
        assert_eq!(ellipf(0.5, NAN), Err("ellipf: Arguments cannot be NAN."));
        // phi = inf: F(inf, m) = inf
        assert_eq!(ellipf(INFINITY, 0.5).unwrap(), INFINITY);
        // phi = -inf: F(-inf, m) = -inf
        assert_eq!(ellipf(NEG_INFINITY, 0.5).unwrap(), NEG_INFINITY);
        // phi > 1/epsilon: F(phi, m) = 2 * phi * K(m) / pi
        assert_eq!(
            ellipf(1e16, 0.4).unwrap(),
            2.0 * 1e16 * ellipk(0.4).unwrap() / PI
        );
        // |phi| -> 0: F(phi, m) = sin(phi)
        assert_eq!(ellipf(1e-100, 0.4).unwrap(), (1e-100).sin());
        assert_eq!(ellipf(-1e-100, 0.4).unwrap(), (-1e-100).sin());
        // m = inf: should return Err
        assert_eq!(
            ellipf(0.5, INFINITY),
            Err("ellipf: m sin²φ must be smaller than one.")
        );
        // m = -inf: F(phi, -inf) = 0.0
        assert_eq!(ellipf(0.5, NEG_INFINITY).unwrap(), 0.0);
    }
}

#[cfg(feature = "test_force_fail")]
crate::test_force_unreachable! {
    assert_eq!(ellipf(0.5, 0.5), Err("ellipf: Unexpected error."));
}
