/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use num_traits::Float;

use crate::{carlson::elliprj_unchecked, crate_util::check, ellipf, ellipk, StrErr};

/// Computes [Heuman Lambda](https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/ellint/heuman_lambda.html).
/// ```text
///            F(φ, 1 - m)   2                     
/// Λ0(φ, m) = ─────────── + ─ ⋅ K(m) ⋅ Z(φ, 1 - m)
///             K(1 - m)     π                     
/// ```
///
/// ## Parameters
/// - phi: amplitude angle (φ). φ ∈ ℝ.
/// - m: elliptic parameter. m ∈ ℝ, m ∈ [0, 1).
///
/// The elliptic modulus (k) is also frequently used instead of the parameter (m), where k² = m.
///
/// ## Domain
/// - Returns error if m < 0 or m ≥ 1.
/// - Returns error if phi is infinite.
///
/// ## Graph
/// ![Heuman Lambda Function](https://github.com/p-sira/ellip/blob/main/figures/heuman_lambda.svg?raw=true)
///
/// [Interactive Plot](https://p-sira.github.io/ellippy/_static/figures/heuman_lambda.html)
///
/// ## Special Cases
/// - Λ0(nπ/2, m) = n where n ∈ ℤ.
/// - Λ0(φ, 0) = sin(φ)
///
/// # Related Functions
/// With mc = 1 - m and Δ² = 1 - mc sin²φ
/// - [heuman_lambda](crate::heuman_lambda)(φ, m) = [ellipf](crate::ellipf)(φ, mc) / [ellipk](crate::ellipk)(mc) + 2/π * [ellipk](crate::ellipk)(m) * [jacobi_zeta](crate::jacobi_zeta)(φ, mc)
/// - [heuman_lambda](crate::heuman_lambda)(φ, m) = 2/π * mc sin(φ)cos(φ)/Δ * [[elliprf](crate::elliprf)(0,mc,1) + m/3Δ² * [elliprj](crate::elliprj)(0,mc,1,1-m/Δ²)]
///
/// # Examples
/// ```
/// use ellip::{heuman_lambda, util::assert_close};
/// use std::f64::consts::FRAC_PI_4;
///
/// assert_close(heuman_lambda(FRAC_PI_4, 0.5).unwrap(), 0.6183811341833665, 1e-15);
/// ```
///
/// # References
/// - Maddock, John, Paul Bristow, Hubert Holin, and Xiaogang Zhang. “Boost Math Library: Special Functions - Elliptic Integrals.” Accessed August 30, 2025. <https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/ellint.html>.
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn heuman_lambda<T: Float>(phi: T, m: T) -> Result<T, StrErr> {
    let ans = heuman_lambda_unchecked(phi, m);
    #[cfg(not(feature = "test_force_fail"))]
    if ans.is_finite() {
        return Ok(ans);
    }
    check!(@nan, heuman_lambda, [phi, m]);
    if m < 0.0 || m >= 1.0 {
        return Err("heuman_lambda: m must satisfy 0.0 ≤ m < 1.0.");
    }
    check!(@inf, heuman_lambda, [phi]);
    Err("heuman_lambda: Unexpected error.")
}

/// Unsafe version of [heuman_lambda].
/// <div class="warning">⚠️ Unstable feature. May subject to changes.</div>
///
/// Undefined behavior with invalid arguments and edge cases.
/// # Known Invalid Cases
/// - m >= 1
/// - m < 0
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn heuman_lambda_unchecked<T: Float>(phi: T, m: T) -> T {
    if m <= 0.0 {
        if m == 0.0 {
            return phi.sin();
        }
        return nan!();
    }

    let n = (phi / pi_2!()).round();
    if (phi - n * pi_2!()).abs() < epsilon!() {
        return n;
    }

    let mc = 1.0 - m;

    let f = ellipf(phi, mc).unwrap_or(nan!());
    let k_m = ellipk(m).unwrap_or(nan!());
    let k_mc = ellipk(mc).unwrap_or(nan!());
    let zeta = jacobi_zeta_unchecked_k(phi, mc, k_mc);

    f / k_mc + k_m * zeta / pi_2!()
}

/// jacobi_zeta_unchecked with K(m) as an argument
///
/// Assume m < 1 and valid K.
#[inline]
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
fn jacobi_zeta_unchecked_k<T: Float>(phi: T, m: T, k: T) -> T {
    let sign = phi.signum();
    let phi = phi.abs();
    let sinp = phi.sin();
    let cosp = phi.cos();

    let mc = 1.0 - m;
    let c2p = cosp * cosp;
    let one_m_ms2p = mc + m * c2p;

    sign * m * sinp * cosp * one_m_ms2p.sqrt() * elliprj_unchecked(0.0, mc, 1.0, one_m_ms2p)
        / (3.0 * k)
}

#[cfg(not(feature = "test_force_fail"))]
#[cfg(all(test, not(feature = "no_std")))]
mod tests {
    use super::*;
    use crate::compare_test_data_wolfram;

    #[test]
    fn test_heuman_lambda() {
        compare_test_data_wolfram!("heuman_lambda_data.csv", heuman_lambda, 2, 3e-15);
    }

    #[test]
    fn test_heuman_lambda_special_cases() {
        use std::f64::{
            consts::{FRAC_PI_2, PI},
            INFINITY, NAN,
        };

        // lambda(nπ/2, m) = n for any integer n
        assert_eq!(heuman_lambda(FRAC_PI_2, 0.5).unwrap(), 1.0);
        assert_eq!(heuman_lambda(PI, 0.5).unwrap(), 2.0);
        assert_eq!(heuman_lambda(3.0 * FRAC_PI_2, 0.5).unwrap(), 3.0);
        // m > 1: should return Err
        assert_eq!(
            heuman_lambda(1.0, 1.5),
            Err("heuman_lambda: m must satisfy 0.0 ≤ m < 1.0.")
        );
        // m = 0: sin(phi)
        assert_eq!(heuman_lambda(1.0, 0.0).unwrap(), 1.0.sin());
        // m < 0: should return Err
        assert_eq!(
            heuman_lambda(1.0, -1.0),
            Err("heuman_lambda: m must satisfy 0.0 ≤ m < 1.0.")
        );
        // NANs: should return Err
        assert_eq!(
            heuman_lambda(NAN, 1.0),
            Err("heuman_lambda: Arguments cannot be NAN.")
        );
        assert_eq!(
            heuman_lambda(1.0, NAN),
            Err("heuman_lambda: Arguments cannot be NAN.")
        );
        // inf: should return Err
        assert_eq!(
            heuman_lambda(INFINITY, 0.5),
            Err("heuman_lambda: phi cannot be infinite.")
        );
        assert_eq!(
            heuman_lambda(-INFINITY, 0.5),
            Err("heuman_lambda: phi cannot be infinite.")
        );
    }
}

#[cfg(feature = "test_force_fail")]
crate::test_force_unreachable! {
    assert_eq!(heuman_lambda(0.5, 0.5), Err("heuman_lambda: Unexpected error."));
}
