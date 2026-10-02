/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Boost Math, see LICENSE in this directory.
 */

// Original header from Boost Math
//  Copyright (c) 2006 Xiaogang Zhang
//  Copyright (c) 2006 John Maddock
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0.

use num_traits::Float;

use crate::{
    carlson::{elliprf_unchecked, elliprj_unchecked},
    crate_util::check,
    ellipe, ellipk, StrErr,
};

/// Computes [complete elliptic integral of the third kind](https://dlmf.nist.gov/19.2.E8).
/// ```text
///              π/2                              
///             ⌠                 dϑ              
/// Π(n, m)  =  ⎮ ──────────────────────────────────
///             ⎮   _____________                
///             ⌡ ╲╱ 1 - m sin²ϑ  ⋅ ( 1 - n sin²ϑ )
///            0              
/// ```
///
/// ## Parameters
/// - n: characteristic, n ∈ ℝ, n ≠ 1.
/// - m: elliptic parameter. m ∈ ℝ, m ≤ 1.
///
/// The elliptic modulus (k) is frequently used instead of the parameter (m), where k² = m.
/// The characteristic (n) is also sometimes expressed in term of α, where α² = n.
///
/// ## Domain
/// - Returns error if n = 1 or m > 1.
/// - Returns the Cauchy principal value if n > 1.
///
/// ## Graph
/// ![Complete Elliptic Integral of the Third Kind](https://github.com/p-sira/ellip/blob/main/figures/ellippi_3d.svg?raw=true)
///
/// [Interactive Plot](https://p-sira.github.io/ellippy/_static/figures/ellippi_3d.html)
///
/// ## Special Cases
/// - Π(0, 0) = π/2
/// - Π(0, m) = K(m)
/// - Π(n, 0) = π/(2 sqrt(1-n)) for n < 1
/// - Π(n, 0) = 0 for n > 1
/// - Π(n, m) = ∞ for n -> 1-
/// - Π(n, 1) = sign(1-n) ∞
/// - Π(∞, m) = Π(-∞, m) = 0
/// - Π(n, -∞) = 0
///
/// # Related Functions
/// - [ellippi](crate::ellippi)(n, m) = n / 3 * [elliprj](crate::elliprj)(0, 1 - m, 1, 1 - n) + [ellipk](crate::ellipk)(m)
/// - [ellippi](crate::ellippi)(n, n) = [ellipe](crate::ellipe)(n) / (1-n) for n < 1
/// - [ellippi](crate::ellippi)(n, m) = [ellipk](crate::ellipk)(m) - ([ellipe](crate::ellipe)(m) / (1-m)) for n -> 1+
///
/// # Examples
/// ```
/// use ellip::{ellippi, util::assert_close};
///
/// assert_close(ellippi(0.5, 0.5).unwrap(), 2.7012877620953506, 1e-15);
/// ```
///
/// # References
/// - Maddock, John, Paul Bristow, Hubert Holin, and Xiaogang Zhang. “Boost Math Library: Special Functions - Elliptic Integrals.” Accessed April 17, 2025. <https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/ellint.html>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn ellippi<T: Float>(n: T, m: T) -> Result<T, StrErr> {
    check!(@nan, ellippi, [n, m]);
    if m > 1.0 {
        return Err("ellippi: m must not be greater than 1.");
    }
    if n == 1.0 {
        return Err("ellippi: n cannot be 1.");
    }
    if m == 1.0 {
        return Ok((1.0 - n).signum() * inf!());
    }
    if n > 1.0 {
        // Use the Cauchy principal value. The n -> 1+ limit is not uniform as
        // m -> 1-, so the direct Carlson form is required at that corner.
        // https://dlmf.nist.gov/19.25.E4
        return Ok(-3.0.recip() * m / n * elliprj_unchecked(0.0, 1.0 - m, 1.0, 1.0 - m / n));
    }

    let ans = ellippi_unchecked(n, m);
    #[cfg(not(feature = "test_force_fail"))]
    if ans.is_finite() {
        return Ok(ans);
    }
    check!(@nan, ellippi, [n, m]);
    if m >= 1.0 {
        if m > 1.0 {
            return Err("ellippi: m must not be greater than 1.");
        }
        // m -> 1-
        let sign = (1.0 - n).signum();
        return Ok(sign * inf!());
    }
    let lim_min = 1e-2 * min_val!();
    if n < lim_min || m < lim_min {
        // n = -inf: Π(-inf, m) = 0
        // m = -inf: Π(n, -inf) = 0
        return Ok(0.0);
    }
    Err("ellippi: Unexpected error.")
}

/// Unsafe version of [ellippi](crate::ellippi).
/// <div class="warning">⚠️ Unstable feature. May subject to changes.</div>
///
/// Undefined behavior with invalid arguments and edge cases.
/// # Known Invalid Cases
/// - n -> 1 or n = 1
/// - n > 1 (p.v. cases)
/// - m -> 1 or m = 1
/// - m > 1
/// - n -> -∞ or m -> -∞
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub fn ellippi_unchecked<T: Float>(n: T, m: T) -> T {
    if n <= 0.0 {
        if m == 1.0 {
            return -inf!();
        }
        if n == 0.0 {
            if m == 0.0 {
                // ellippi(0,0)
                return pi_2!();
            }
            // ellippi(0,m)
            return ellipk(m).unwrap_or(nan!());
        }
        // n < 0: ellippi(n,m)
        // When n < 0 and m < 0, standard A&S 17.7.17 suffers from catastrophic cancellation in m - n.
        // Substituting Π((m-n)/(1-n), m) = K(m) + (m-n)/(3(1-n)) R_J into A&S 17.7.17:
        //   Π(n, m) = -n(1-m)/((1-n)(m-n)) [K(m) + (m-n)/(3(1-n)) R_J] + m/(m-n) K(m)
        // The (m - n) denominator cancels out algebraically:
        //   [m/(m-n) - n(1-m)/((1-n)(m-n))] K(m) = 1/(1-n) K(m)
        // Leaving the cancellation-free, purely additive form:
        //   Π(n, m) = 1/(1-n) K(m) - n(1-m)/(3(1-n)²) R_J(0, 1-m, 1, (1-m)/(1-n))
        // Valid for all n < 0 and m < 0, eliminating division by (m - n) and near-diagonal cancellation.
        if m < 0.0 {
            if m == n {
                return ellipe(m).unwrap_or(nan!()) / (1.0 - m);
            }
            let one_minus_n = 1.0 - n;
            let one_minus_m = 1.0 - m;
            let p = one_minus_m / one_minus_n;
            let km = ellipk(m).unwrap_or(nan!());
            let rj = elliprj_unchecked(0.0, one_minus_m, 1.0, p);
            return km / one_minus_n + (-n) * one_minus_m / (3.0 * one_minus_n * one_minus_n) * rj;
        }

        // Apply A&S 17.7.17 for n < 0 and m >= 0
        let nn = (m - n) / (1.0 - n);
        let nm1 = (1.0 - m) / (1.0 - n);

        let mut result = ellippi_vc(nn, m, nm1);
        // Split calculations to avoid overflow/underflow
        result = result * -n / (1.0 - n);
        result = result * (1.0 - m) / (m - n);
        result = result + ellipk(m).unwrap_or(nan!()) * m / (m - n);
        return result;
    }

    // https://dlmf.nist.gov/19.6.E1
    if m == n {
        let mc = 1.0 - m;
        return 1.0 / mc * ellipe(m).unwrap_or(nan!());
    }

    // Compute vc = 1-n without cancellation errors
    let vc = 1.0 - n;
    ellippi_vc(n, m, vc)
}

#[inline]
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn ellippi_vc<T: Float>(n: T, m: T, vc: T) -> T {
    let x = 0.0;
    let y = 1.0 - m;
    let z = 1.0;
    let p = vc;

    elliprf_unchecked(x, y, z) + n * elliprj_unchecked(x, y, z, p) / 3.0
}

#[cfg(not(feature = "test_force_fail"))]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{compare_test_data_boost, compare_test_data_wolfram};

    #[test]
    fn test_ellippi() {
        compare_test_data_boost!("ellippi2_data_f64.txt", ellippi, 2, 5e-16);
    }

    #[test]
    fn test_ellippi_wolfram() {
        compare_test_data_wolfram!("ellippi_data.csv", ellippi, 2, 5e-15);
    }

    #[test]
    fn test_ellippi_wolfram_neg() {
        // A16 removes the stale dedicated dataset test because these cases now live here:
        // https://github.com/p-sira/ellip/pull/130
        // Includes the diagonal/near-diagonal cases generated by ellippi.wls.
        compare_test_data_wolfram!("ellippi_neg.csv", ellippi, 2, 1.3e-14);
    }

    #[test]
    fn test_ellippi_wolfram_pv() {
        compare_test_data_wolfram!("ellippi_pv.csv", ellippi, 2, 1e-14);
    }

    #[test]
    fn test_ellippi_negative_diagonal() {
        // Pi(n, n) = E(n) / (1 - n) for n < 1. Previously errored for n = m < 0 because
        // the A&S 17.7.17 branch divides by m - n. Reference values from mpmath.
        use crate::util::assert_close;
        assert_close(ellippi(-0.5, -0.5).unwrap(), 1.1678475171298786, 1e-14);
        assert_close(ellippi(-2.0, -2.0).unwrap(), 0.72814604758206706, 1e-14);
        assert_close(ellippi(-10.0, -10.0).unwrap(), 0.33083073076525165, 1e-14);
        assert_close(ellippi(-100.0, -100.0).unwrap(), 0.10108179128529279, 1e-14);
        assert_close(
            ellippi(-1000.0, -1000.0).unwrap(),
            0.031675528525186074,
            1e-14,
        );
        // Consistency with the closed form and with the (already handled) n = m > 0 diagonal.
        assert_close(
            ellippi(-0.5, -0.5).unwrap(),
            ellipe(-0.5).unwrap() / 1.5,
            1e-14,
        );
        assert_close(
            ellippi(0.5, 0.5).unwrap(),
            ellipe(0.5).unwrap() / 0.5,
            1e-14,
        );
    }

    #[test]
    fn test_ellippi_special_cases() {
        use std::f64::{
            consts::{FRAC_PI_2, PI},
            EPSILON, INFINITY, NAN, NEG_INFINITY,
        };
        // m > 1: should return Err
        assert_eq!(
            ellippi(0.5, 1.1),
            Err("ellippi: m must not be greater than 1.")
        );
        // n == 1: should return Err
        assert_eq!(ellippi(1.0, 0.5), Err("ellippi: n cannot be 1."));
        // n = 0: Π(0, m) = K(m)
        assert_eq!(ellippi(0.0, 0.5).unwrap(), ellipk(0.5).unwrap());
        // m = 0, n < 1: Π(n, 0) = pi/(2 sqrt(1-n))
        assert_eq!(ellippi(0.5, 0.0).unwrap(), PI / (2.0 * 0.5.sqrt()));
        // m = 0, n > 1: Π(n, 0) = 0
        assert_eq!(ellippi(2.0, 0.0).unwrap(), 0.0);
        // m = 1: Π(n, 1) = sign(1-n) inf
        assert_eq!(ellippi(2.0, 1.0).unwrap(), NEG_INFINITY);
        assert_eq!(ellippi(0.5, 1.0).unwrap(), INFINITY);
        assert_eq!(ellippi(-2.0, 1.0).unwrap(), INFINITY);
        // The last representable n below 1 is still finite (100-digit reference).
        crate::util::assert_close(
            ellippi(1.0 - 0.5 * EPSILON, 0.5).unwrap(),
            210828713.28594348,
            2e-15,
        );
        // The nearest representable n > 1 remains finite.
        crate::util::assert_close(
            ellippi(1.0 + EPSILON, 0.5).unwrap(),
            -0.8472130847939789,
            2e-15,
        );
        // Π(0, 0) = pi/2
        assert_eq!(ellippi(0.0, 0.0).unwrap(), FRAC_PI_2);
        // Π(n, n) = E(n) / (1-n) for n < 1
        assert_eq!(ellippi(0.5, 0.5).unwrap(), ellipe(0.5).unwrap() / 0.5);
        // n = inf: Π(inf, m) = 0
        assert_eq!(ellippi(INFINITY, 0.5).unwrap(), 0.0);
        // n = -inf: Π(-inf, m) = 0
        assert_eq!(ellippi(NEG_INFINITY, 0.5).unwrap(), 0.0);
        // m = -inf: Π(n, -inf) = 0
        assert_eq!(ellippi(0.5, NEG_INFINITY).unwrap(), 0.0);
        // n = nan or m = nan: should return Err
        assert_eq!(ellippi(NAN, 0.5), Err("ellippi: Arguments cannot be NAN."));
        assert_eq!(ellippi(0.5, NAN), Err("ellippi: Arguments cannot be NAN."));
        // m = inf: should return Err
        assert_eq!(
            ellippi(0.5, INFINITY),
            Err("ellippi: m must not be greater than 1.")
        );
    }

    // Regression for audit finding A5: https://github.com/p-sira/ellip/pull/119
    #[test]
    fn test_ellippi_invalid_input_does_not_abort() {
        assert!(ellippi(f64::NAN, 0.5).is_err());
        assert!(ellippi(2.0, f64::NAN).is_err());
        if std::env::var("ELLIP_A05_PI_CHILD").is_ok() {
            assert!(ellippi(2.0, 2.0).is_err());
            return;
        }
        let status = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "legendre::ellippi::tests::test_ellippi_invalid_input_does_not_abort",
            ])
            .env("ELLIP_A05_PI_CHILD", "1")
            .status()
            .unwrap();
        assert!(
            status.success(),
            "invalid Pi input aborted or failed: {status}"
        );
    }

    // Regression for audit finding A4: https://github.com/p-sira/ellip/pull/118
    #[test]
    fn test_finite_below_singular_boundary() {
        let u = f64::from_bits(1.0f64.to_bits() - 1);
        crate::assert_close!(ellippi(u, 0.0).unwrap(), 149078413.4323951, 2e-15);
        crate::assert_close!(ellippi(0.0, u).unwrap(), 19.75469464595844, 2e-15);
        assert!(ellippi(1.0, 0.0).is_err());
        assert_eq!(ellippi(0.0, 1.0).unwrap(), f64::INFINITY);
        let u = f32::from_bits(1.0f32.to_bits() - 1);
        crate::assert_close!(
            ellippi(u, 0.0).unwrap() as f64,
            (std::f32::consts::FRAC_PI_2 / (1.0 - u).sqrt()) as f64,
            4e-7
        );
        crate::assert_close!(
            ellippi(0.0, u).unwrap() as f64,
            ellipk(u).unwrap() as f64,
            4e-7
        );
    }

    // Regression for branch-limit margin analysis: https://github.com/p-sira/ellip/pull/132
    #[test]
    fn test_principal_value_at_simultaneous_branch_limits() {
        // Exact-binary, 80-digit Wolfram references. This corner used the nonuniform
        // n -> 1+ limit and lost most significant digits when m also approached 1.
        let n = f64::from_bits(1.0f64.to_bits() + 1);
        let m = f64::from_bits(1.0f64.to_bits() - 1);
        crate::assert_close!(ellippi(n, m).unwrap(), -4_214_834_719_445_440.5, 3e-15);

        let n = f32::from_bits(1.0f32.to_bits() + 1);
        let m = f32::from_bits(1.0f32.to_bits() - 1);
        crate::assert_close!(ellippi(n, m).unwrap() as f64, -7_850_737.0, 4e-7);
    }
}

#[cfg(feature = "test_force_fail")]
crate::test_force_unreachable! {
    assert_eq!(ellippi(0.5, 0.5), Err("ellippi: Unexpected error."));
}
