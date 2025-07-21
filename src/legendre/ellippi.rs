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

use crate::{crate_util::check_nan, ellipe, ellipk, elliprf, elliprj, StrErr};

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
/// - m: elliptic parameter. m ∈ ℝ, m < 1.
///
/// The elliptic modulus (k) is frequently used instead of the parameter (m), where k² = m.
/// The characteristic (n) is also sometimes expressed in term of α, where α² = n.
///
/// ## Domain
/// - Returns error if n = 1 or m > 1.
/// - Returns the Cauchy principal value if n > 1.
///
/// ## Graph
/// ![Complete Elliptic Integral of the Third Kind](https://github.com/p-sira/ellip/blob/main/figures/ellippi_plot.svg?raw=true)
///
/// [Interactive Plot](https://github.com/p-sira/ellip/blob/main/figures/ellippi_plot.html)
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
    check_nan!(ellippi, [n, m]);

    if m > 1.0 {
        return Err("ellippi: m must be less than 1.");
    }

    if n == 1.0 {
        return Err("ellippi: n cannot be 1.");
    }

    // n = -inf: Π(-inf, m) = 0
    // m = -inf: Π(n, -inf) = 0
    if n == neg_inf!() || m == neg_inf!() {
        return Ok(0.0);
    }

    // m -> 1-
    if 1.0 - m <= epsilon!() {
        let sign = (1.0 - n).signum();
        return Ok(sign * inf!());
    }

    if n > 1.0 {
        // n -> 1+
        // https://dlmf.nist.gov/19.6.E6
        if n - 1.0 <= epsilon!() {
            return Ok(ellipk(m)? - ellipe(m)? / (1.0 - m));
        }

        // Use Cauchy principal value
        // https://dlmf.nist.gov/19.25.E4
        return Ok(-1.0 / 3.0 * m / n * elliprj(0.0, 1.0 - m, 1.0, 1.0 - m / n)?);
    }

    // n < 1 and n -> 1-
    if 1.0 - n <= epsilon!() {
        return Ok(inf!());
    }

    if n == 0.0 {
        if m == 0.0 {
            return Ok(pi_2!());
        }
        return ellipk(m);
    }

    // https://dlmf.nist.gov/19.6.E1
    if m == n {
        let mc = 1.0 - m;
        return Ok(1.0 / mc * ellipe(m)?);
    }

    if n < 0.0 {
        // Apply A&S 17.7.17
        let nn = (m - n) / (1.0 - n);
        let nm1 = (1.0 - m) / (1.0 - n);

        let mut result = ellippi_vc(nn, m, nm1)?;
        // Split calculations to avoid overflow/underflow
        result = result * -n / (1.0 - n);
        result = result * (1.0 - m) / (m - n);
        result = result + ellipk(m)? * m / (m - n);
        return Ok(result);
    }

    // Compute vc = 1-n without cancellation errors
    let vc = 1.0 - n;
    ellippi_vc(n, m, vc)
}

#[inline]
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn ellippi_vc<T: Float>(n: T, m: T, vc: T) -> Result<T, StrErr> {
    let x = 0.0;
    let y = 1.0 - m;
    let z = 1.0;
    let p = vc;

    Ok(elliprf(x, y, z)? + n * elliprj(x, y, z, p)? / 3.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compare_test_data_boost;

    fn ellippi_k(inp: &[f64]) -> f64 {
        ellippi(inp[0], inp[1] * inp[1]).unwrap()
    }

    #[test]
    fn test_ellippi() {
        compare_test_data_boost!("ellippi2_data_f64.txt", ellippi_k, 4.4e-16);
    }

    #[test]
    fn test_ellippi_err() {
        // m > 1
        assert!(ellippi(0.5, 1.1).is_err());
        // n == 1
        assert!(ellippi(1.0, 0.5).is_err());
    }

    #[test]
    fn test_ellippi_special_cases() {
        use std::f64::{
            consts::{FRAC_PI_2, PI},
            EPSILON, INFINITY, NAN, NEG_INFINITY,
        };
        // m > 1: should return Err
        assert!(ellippi(0.5, 1.1).is_err());
        // n == 1: should return Err
        assert!(ellippi(1.0, 0.5).is_err());
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
        // n -> 1-: Π(n, m) = inf
        assert_eq!(ellippi(1.0 - EPSILON, 0.5).unwrap(), INFINITY);
        // n -> 1+: Π(n, m) = K(m) - E(m) / (1-m)
        assert_eq!(
            ellippi(1.0 + EPSILON, 0.5).unwrap(),
            ellipk(0.5).unwrap() - ellipe(0.5).unwrap() / 0.5
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
        assert!(ellippi(NAN, 0.5).is_err());
        assert!(ellippi(0.5, NAN).is_err());
        // m = inf: should return Err
        assert!(ellippi(0.5, INFINITY).is_err());
    }
}
