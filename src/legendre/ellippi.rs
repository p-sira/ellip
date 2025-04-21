/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Boost Math, see LICENSE in this directory.
 */

// Original header from Boost Math
//  Copyright (c) 2006 Xiaogang Zhang
//  Copyright (c) 2006 John Maddock
//  Use, modification and distribution are subject to the
//  Boost Software License, Version one!(). (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  History:
//  XZ wrote the original of this file as part of the Google
//  Summer of Code 2006.  JM modified it to fit into the
//  Boost.Math conceptual framework better, and to correctly
//  handle the various corner cases.

use num_traits::Float;

use crate::{ellipe, ellipk, elliprf, elliprj};

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
/// - Returns the principal value if n > 1.
///
/// ## Graph
/// ![Complete Elliptic Integral of the Third Kind](https://github.com/p-sira/ellip/blob/main/figures/ellippi_plot.svg?raw=true)
///
/// [Interactive Plot](https://github.com/p-sira/ellip/blob/main/figures/ellippi_plot.html)
///
/// # Related Functions
/// - [ellippi](crate::ellippi)(n, m) = n / 3 * [elliprj](crate::elliprj)(0, 1 - m, 1, 1 - n) + [ellipk](crate::ellipk)(m)
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
///
pub fn ellippi<T: Float>(n: T, m: T) -> Result<T, &'static str> {
    if m > one!() {
        return Err("ellippi: m must be less than 1.");
    }

    if n == one!() {
        return Err("ellippi: n cannot be 1.");
    }

    // m -> 1-
    if one!() - m <= epsilon!() {
        let sign = (one!() - n).signum();
        return Ok(sign * inf!());
    }

    if n > one!() {
        // n -> 1+
        // https://dlmf.nist.gov/19.6.E6
        if n - one!() <= epsilon!() {
            return Ok(ellipk(m)? - ellipe(m)? / (one!() - m));
        }

        // Use Cauchy principal value
        // https://dlmf.nist.gov/19.25.E4
        return Ok(-third!() * m / n * elliprj(zero!(), one!() - m, one!(), one!() - m / n)?);
    }

    // n < 1 and n -> 1-
    if one!() - n <= epsilon!() {
        return Ok(inf!());
    }

    if n == zero!() {
        if m == zero!() {
            return Ok(pi_2!());
        }
        return ellipk(m);
    }

    if n < zero!() {
        // When m < 0, n < 0 and m == n, Boost implementation cancels out, resulting in inf.
        // https://dlmf.nist.gov/19.6.E13 with phi = π/2
        if m == n {
            let mc = one!() - m;
            return Ok(one!() / mc * ellipe(m)?);
        }

        // Apply A&S 17.7.17
        let nn = (m - n) / (one!() - n);
        let nm1 = (one!() - m) / (one!() - n);

        let mut result = ellippi_vc(nn, m, nm1)?;
        // Split calculations to avoid overflow/underflow
        result = result * -n / (one!() - n);
        result = result * (one!() - m) / (m - n);
        result = result + ellipk(m)? * m / (m - n);
        return Ok(result);
    }

    // Compute vc = 1-n without cancellation errors
    let vc = one!() - n;
    ellippi_vc(n, m, vc)
}

#[inline]
pub fn ellippi_vc<T: Float>(n: T, m: T, vc: T) -> Result<T, &'static str> {
    let x = zero!();
    let y = one!() - m;
    let z = one!();
    let p = vc;

    Ok(elliprf(x, y, z)? + n * elliprj(x, y, z, p)? / three!())
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
        compare_test_data_boost!(
            "./tests/data/boost/ellippi2_data_f64.txt",
            ellippi_k,
            4.4e-16
        );
    }
}
