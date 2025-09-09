/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Boost Math and SciPy, see LICENSE in this directory.
 */

// Original header from Boost Math
//  Copyright (c) 2006 Xiaogang Zhang
//  Copyright (c) 2006 John Maddock
//  Copyright (c) 2024 Matt Borland
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0.

use num_traits::Float;

use crate::{
    carlson::{elliprd_unchecked, elliprf_unchecked},
    crate_util::check,
    ellipe, StrErr,
};

/// Computes [incomplete elliptic integral of the second kind](https://dlmf.nist.gov/19.2.E5).
/// ```text
///              φ
///             ⌠   _____________
/// E(φ, m)  =  │ \╱ 1 - m sin²θ  dθ
///             ⌡
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
/// ![Incomplete Elliptic Integral of the Second Kind](https://github.com/p-sira/ellip/blob/main/figures/ellipeinc_plot.svg?raw=true)
///
/// [Interactive Plot](https://github.com/p-sira/ellip/blob/main/figures/ellipeinc_plot.html)
///
/// ![Incomplete Elliptic Integral of the Second Kind (3D Plot)](https://github.com/p-sira/ellip/blob/main/figures/ellipeinc_plot_3d.png?raw=true)
///
/// [Interactive 3D Plot](https://github.com/p-sira/ellip/blob/main/figures/ellipeinc_plot_3d.html)
///
/// ## Special Cases
/// - E(0, m) = 0
/// - E(φ, 0) = φ
/// - E(π/2, m) = E(m)
/// - E(φ, 1) = sin(φ)
/// - E(φ, -∞) = ∞
/// - E(∞, m) = ∞
/// - E(-∞, m) = -∞
///
/// # Related Functions
/// With c = csc²φ,
/// - [ellipeinc](crate::ellipeinc)(φ, m) = [elliprf](crate::elliprf)(c - 1, c - m, c) - m / 3 * [elliprd](crate::elliprd)(c - 1, c - m, c)
///
/// # Examples
/// ```
/// use ellip::{ellipeinc, util::assert_close};
/// use std::f64::consts::FRAC_PI_4;
///
/// assert_close(ellipeinc(FRAC_PI_4, 0.5).unwrap(), 0.7481865041776612, 1e-15);
/// ```
///
/// # References
/// - Maddock, John, Paul Bristow, Hubert Holin, and Xiaogang Zhang. “Boost Math Library: Special Functions - Elliptic Integrals.” Accessed April 17, 2025. <https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/ellint.html>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
/// - The MathWorks, Inc. “ellipticE.” Accessed April 21, 2025. <https://www.mathworks.com/help/symbolic/sym.elliptice.html>.
pub fn ellipeinc<T: Float>(phi: T, m: T) -> Result<T, StrErr> {
    let ans = ellipeinc_unchecked(phi, m)?;
    // Infinites are expected to be handled properly by ellipeinc_unchecked
    if !ans.is_nan() {
        return Ok(ans);
    }
    #[cfg(feature = "test_force_fail")]
    let ans = nan!();
    check!(@nan, ellipeinc, [phi, m]);
    Err("ellipeinc: Unexpected error.")
}

/// Unsafe version of [ellipeinc].
/// <div class="warning">⚠️ Unstable feature. May subject to changes.</div>
///
/// Undefined behavior with invalid arguments and edge cases.
/// # Known Invalid Cases
/// - NAN arguments
#[inline]
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn ellipeinc_unchecked<T: Float>(phi: T, m: T) -> Result<T, StrErr> {
    if phi == 0.0 {
        return Ok(0.0);
    }

    if m == 0.0 {
        return Ok(phi);
    }

    let (phi, invert) = if phi < 0.0 { (-phi, -1.0) } else { (phi, 1.0) };

    if phi > 1.0 / epsilon!() {
        if phi >= max_val!() {
            return Ok(invert * inf!());
        }
        return Ok(invert * 2.0 * phi * ellipe_wrapper(m)? / pi!());
    }

    // m -> -inf, sqrt(1+m sin²θ) -> m sinθ
    // E(phi, m) -> sqrt(-m) int(sin(theta)) d(theta), theta=[0,phi]
    // E(phi, m) -> sqrt(-m) (1-cos(phi))
    if m <= -max_val!() {
        return Ok((1.0 - phi.cos()) * (-m).sqrt());
    }

    if m == 1.0 {
        // For k = 1 ellipse actually turns to a line and every pi/2 in phi is exactly 1 in arc length
        // Periodicity though is in pi, curve follows sin(pi) for 0 <= phi <= pi/2 and then
        // 2 - sin(pi- phi) = 2 + sin(phi - pi) for pi/2 <= phi <= pi, so general form is:
        //
        // 2n + sin(phi - n * pi) ; |phi - n * pi| <= pi / 2
        let mm = (phi / pi!()).round();
        let remains = phi - mm * pi!();
        return Ok(invert * (2.0 * mm + remains.sin()));
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

    let mut result = if m > 0.0 && rphi.powi(3) * m / 6.0 < epsilon!() * rphi.abs() {
        // See http://functions.wolfram.com/EllipticIntegrals/EllipticE2/06/01/03/0001/
        s * rphi
    } else {
        let s2p = rphi.sin() * rphi.sin();
        if m * s2p >= 1.0 {
            return Err("ellipeinc: m sin²φ must be smaller than one.");
        }
        let c2p = rphi.cos() * rphi.cos();
        let c = 1.0 / s2p;
        let cm1 = c2p / s2p;
        s * ((1.0 - m) * elliprf_unchecked(cm1, c - m, c)
            + m * (1.0 - m) * elliprd_unchecked(cm1, c, c - m) / 3.0
            + m * (cm1 / (c * (c - m))).sqrt())
    };

    if mm != 0.0 {
        result = result + mm * ellipe_wrapper(m)?;
    }

    Ok(invert * result)
}

#[inline]
fn ellipe_wrapper<T: Float>(m: T) -> Result<T, StrErr> {
    check!(@nan, ellipeinc, [m]);
    ellipe(m)
}

#[cfg(not(feature = "test_force_fail"))]
#[cfg(all(test, not(feature = "no_std")))]
mod tests {
    use crate::compare_test_data_boost;

    use super::*;

    #[test]
    fn test_ellipeinc() {
        compare_test_data_boost!("ellipeinc_data.txt", ellipeinc, 2, 5e-16);
    }

    #[test]
    fn test_ellipeinc_special_cases() {
        use std::f64::{
            consts::{FRAC_PI_2, PI},
            INFINITY, MAX, MIN_POSITIVE, NAN, NEG_INFINITY,
        };
        // phi = pi/2, m=1: E(pi/2, 1) = 1
        assert_eq!(ellipeinc(FRAC_PI_2, 1.0).unwrap(), 1.0);
        // m = 1: E(phi, 1) = sin(phi)
        assert_eq!(ellipeinc(0.4, 1.0).unwrap(), 0.4.sin());
        // m * sin^2(phi) >= 1: should return Err
        assert_eq!(
            ellipeinc(FRAC_PI_2, 2.0),
            Err("ellipeinc: m sin²φ must be smaller than one.")
        );
        // phi = 0: E(0, m) = 0
        assert_eq!(ellipeinc(0.0, 0.5).unwrap(), 0.0);
        // phi -> 0, m > 0: E(phi, m) = phi
        assert_eq!(ellipeinc(MIN_POSITIVE, 0.5).unwrap(), MIN_POSITIVE);
        // phi = pi/2, m = 0: E(pi/2, 0) = pi/2
        assert_eq!(ellipeinc(FRAC_PI_2, 0.0).unwrap(), FRAC_PI_2);
        // m < 0: should be valid
        assert!(ellipeinc(FRAC_PI_2, -1.0).unwrap().is_finite());
        // phi = nan or m = nan: should return Err
        assert_eq!(
            ellipeinc(NAN, 0.5),
            Err("ellipeinc: Arguments cannot be NAN.")
        );
        assert_eq!(
            ellipeinc(0.5, NAN),
            Err("ellipeinc: Arguments cannot be NAN.")
        );
        // phi > 1/epsilon: E(phi, m) = 2 * phi * E(m) / pi
        assert_eq!(
            ellipeinc(1e16, 0.5).unwrap(),
            2.0 * 1e16 * ellipe(0.5).unwrap() / PI
        );
        // phi = inf: E(inf, m) = inf
        assert_eq!(ellipeinc(INFINITY, 0.5).unwrap(), INFINITY);
        // phi = -inf: E(-inf, m) = -inf
        assert_eq!(ellipeinc(NEG_INFINITY, 0.5).unwrap(), NEG_INFINITY);
        // m = inf: should return Err
        assert_eq!(
            ellipeinc(0.5, INFINITY),
            Err("ellipeinc: m sin²φ must be smaller than one.")
        );
        // m -> -inf: E(phi, m) = (1-cos(phi)) sqrt(-m)
        assert_eq!(
            ellipeinc(0.5, -MAX).unwrap(),
            (1.0 - 0.5.cos()) * MAX.sqrt()
        );
        // m = -inf: E(phi, -inf) = inf
        assert_eq!(ellipeinc(0.5, NEG_INFINITY).unwrap(), INFINITY);
    }
}

#[cfg(feature = "test_force_fail")]
crate::test_force_unreachable! {
    assert_eq!(ellipeinc(0.5, 0.2), Err("ellipeinc: Unexpected error."));
}
