/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use num_traits::Float;

use crate::{carlson::elliprj_unchecked, crate_util::check, ellipk, StrErr};

/// Computes [Jacobi Zeta](https://dlmf.nist.gov/22.16.E33).
/// ```text
/// Z(φ, m) = E(φ, m) - E(m) F(φ, m) / K(m)
/// ```
///
/// ## Parameters
/// - phi: amplitude angle (φ). φ ∈ ℝ.
/// - m: elliptic parameter. m ∈ ℝ, m ≤ 1.
///
/// The elliptic modulus (k) is also frequently used instead of the parameter (m), where k² = m.
///
/// ## Domain
/// - Returns error if m > 1.
/// - Returns error if phi or m is infinite.
///
/// ## Graph
/// ![Jacobi Zeta Function](https://github.com/p-sira/ellip/blob/main/figures/jacobi_zeta.svg?raw=true)
///
/// [Interactive Plot](https://p-sira.github.io/ellippy/_static/figures/jacobi_zeta.html)
///
/// ## Special Cases
/// - Z(0, m) = 0
/// - Z(φ, 0) = 0
/// - Z(φ, 1) = sin(φ) sign(cos(φ)) for φ ≠ nπ/2
/// - Z(nπ/2, m) = 0 where n ∈ ℤ.
/// - Z(φ, m) = -Z(-φ, m)
///
/// # Related Functions
/// - [jacobi_zeta](crate::jacobi_zeta)(φ, m) = [ellipeinc](crate::ellipeinc)(φ, m) - [ellipe](crate::ellipe)(m) [ellipf](crate::ellipf)(φ, m) / [ellipk](crate::ellipk)(m)
/// - [jacobi_zeta](crate::jacobi_zeta)(φ, m) = m sin(φ) cos(φ) √(1 - m sin²φ) [elliprj](crate::elliprj)(0, kc², 1, 1 - m sin²φ) / (3 [ellipk](crate::ellipk)(m))
///
/// # Examples
/// ```
/// use ellip::{jacobi_zeta, util::assert_close};
/// use std::f64::consts::FRAC_PI_4;
///
/// assert_close(jacobi_zeta(FRAC_PI_4, 0.5).unwrap(), 0.146454543836188, 1e-15);
/// ```
///
/// # References
/// - Maddock, John, Paul Bristow, Hubert Holin, and Xiaogang Zhang. “Boost Math Library: Special Functions - Elliptic Integrals.” Accessed August 30, 2025. <https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/ellint.html>.
/// - Reinhardt, W. P., and P. L. Walker. “DLMF: Chapter 22 Jacobian Elliptic Functions.” Accessed August 31, 2025. <https://dlmf.nist.gov/22>.
/// - Weisstein, Eric W. “Jacobi Zeta Function.” Wolfram Research, Inc. Accessed August 31, 2025. <https://mathworld.wolfram.com/JacobiZetaFunction.html>.
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn jacobi_zeta<T: Float>(phi: T, m: T) -> Result<T, StrErr> {
    let ans = jacobi_zeta_unchecked(phi, m)?;
    if ans.is_finite() {
        return Ok(ans);
    }
    check!(@nan, jacobi_zeta, [phi, m]);
    check!(@inf, jacobi_zeta, [phi, m]);
    Err("jacobi_zeta: Unexpected error.")
}

/// Unsafe version of [jacobi_zeta].
/// <div class="warning">⚠️ Unstable feature. May subject to changes.</div>
///
/// Undefined behavior with invalid arguments and edge cases.
/// # Known Invalid Cases
/// - |φ| = ∞
/// - m = -∞
#[inline]
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn jacobi_zeta_unchecked<T: Float>(phi: T, m: T) -> Result<T, StrErr> {
    let sign = phi.signum();
    let phi = phi.abs();

    let sinp = phi.sin();
    let cosp = phi.cos();

    let nphi = (phi / (pi_2!())).round();
    Ok(if (phi - nphi * pi_2!()).abs() < epsilon!().sqrt() {
        // Z(nπ/2, m) when n is an integer.
        0.0
    } else if m >= 1.0 {
        if m != 1.0 {
            return Err("jacobi_zeta: m must not be greater than 1.");
        }
        sign * sinp * cosp.signum()
    } else {
        let mc = 1.0 - m;
        let c2p = cosp * cosp;
        let one_m_ms2p = mc + m * c2p;

        sign * m * sinp * cosp * one_m_ms2p.sqrt() * elliprj_unchecked(0.0, mc, 1.0, one_m_ms2p)
            / (3.0 * ellipk(m).unwrap_or(nan!()))
    })
}

#[cfg(not(feature = "test_force_fail"))]
#[cfg(all(test, not(feature = "no_std")))]
mod tests {
    use super::*;
    use crate::compare_test_data_wolfram;

    #[test]
    fn test_jacobi_zeta() {
        compare_test_data_wolfram!("jacobi_zeta_data.csv", jacobi_zeta, 2, 3e-15);
    }

    #[test]
    fn test_jacobi_neg() {
        compare_test_data_wolfram!("jacobi_zeta_neg.csv", jacobi_zeta, 2, 3e-15);
    }

    #[test]
    fn test_jacobi_zeta_special_cases() {
        use std::f64::{
            consts::{FRAC_PI_2, FRAC_PI_3, FRAC_PI_4, FRAC_PI_6, PI},
            INFINITY, NAN,
        };

        // Z(0, m) = 0
        assert_eq!(jacobi_zeta(0.0, 0.5).unwrap(), 0.0);
        assert_eq!(jacobi_zeta(0.0, 0.0).unwrap(), 0.0);
        assert_eq!(jacobi_zeta(0.0, 0.9).unwrap(), 0.0);
        // Z(φ, 0) = 0
        assert_eq!(jacobi_zeta(FRAC_PI_4, 0.0).unwrap(), 0.0);
        assert_eq!(jacobi_zeta(FRAC_PI_2, 0.0).unwrap(), 0.0);
        assert_eq!(jacobi_zeta(PI, 0.0).unwrap(), 0.0);
        // Z(nπ/2, m) = 0 for any integer n
        assert_eq!(jacobi_zeta(FRAC_PI_2, 0.5).unwrap(), 0.0);
        assert_eq!(jacobi_zeta(PI, 0.5).unwrap(), 0.0);
        assert_eq!(jacobi_zeta(3.0 * FRAC_PI_2, 0.5).unwrap(), 0.0);
        // Z(φ, 1) = sin(φ) sign(cos(φ)) for φ ≠ nπ/2
        assert_eq!(jacobi_zeta(FRAC_PI_3, 1.0).unwrap(), FRAC_PI_3.sin());
        assert_eq!(jacobi_zeta(FRAC_PI_4, 1.0).unwrap(), FRAC_PI_4.sin());
        assert_eq!(jacobi_zeta(FRAC_PI_6, 1.0).unwrap(), FRAC_PI_6.sin());
        // m > 1: should return Err
        assert_eq!(
            jacobi_zeta(1.0, 1.5),
            Err("jacobi_zeta: m must not be greater than 1.")
        );
        // NANs: should return Err
        assert_eq!(
            jacobi_zeta(NAN, 1.0),
            Err("jacobi_zeta: Arguments cannot be NAN.")
        );
        assert_eq!(
            jacobi_zeta(1.0, NAN),
            Err("jacobi_zeta: Arguments cannot be NAN.")
        );
        // inf: should return Err
        assert_eq!(
            jacobi_zeta(INFINITY, 1.0),
            Err("jacobi_zeta: phi cannot be infinite.")
        );
        assert_eq!(
            jacobi_zeta(-INFINITY, 1.0),
            Err("jacobi_zeta: phi cannot be infinite.")
        );
        assert_eq!(
            jacobi_zeta(1.0, -INFINITY),
            Err("jacobi_zeta: m cannot be infinite.")
        );
    }
}

#[cfg(feature = "test_force_fail")]
crate::test_force_unreachable! {
    assert_eq!(jacobi_zeta(0.5, 0.5), Err("jacobi_zeta: Unexpected error."));
}
