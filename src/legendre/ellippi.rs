/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Boost Math, see LICENSE in this directory.
 */

// Original header from Boost Math
//  Copyright (c) 2006 Xiaogang Zhang
//  Copyright (c) 2006 John Maddock
//  Use, modification and distribution are subject to the
//  Boost Software License, Version T::one(). (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  History:
//  XZ wrote the original of this file as part of the Google
//  Summer of Code 2006.  JM modified it to fit into the
//  Boost.Math conceptual framework better, and to correctly
//  handle the various corner cases.

use std::f64::consts::FRAC_PI_2;

use num_traits::Float;

use crate::unchecked::{_ellipk, _elliprf, _elliprj};

/// Compute [complete elliptic integral of the third kind](https://dlmf.nist.gov/19.2.E8).
/// ```text
///              π/2                              
///             ⌠                 dϑ              
/// Π(n, m)  =  ⎮ ──────────────────────────────────
///             ⎮   _____________                
///             ⌡ ╲╱ 1 - m sin²ϑ  ⋅ ( 1 - n sin²ϑ )
///            0              
/// where m < 1, n < 1               
/// ```
///
/// Note that some mathematical references use the parameter k and α for the function,
/// where k² = m, α² = n.
pub fn ellippi<T: Float>(n: T, m: T) -> Result<T, &'static str> {
    if m >= T::one() {
        return Err("ellippi: m must be less than 1.");
    }
    if n > T::one() {
        return Err("ellippi: n must be less than 1.");
    }

    if n == T::zero() {
        if m == T::zero() {
            return Ok(T::from(FRAC_PI_2).unwrap());
        }
        return Ok(_ellipk(m));
    }

    if n < T::zero() {
        // Apply A&S 17.7.17
        let nn = (m - n) / (T::one() - n);
        let nm1 = (T::one() - m) / (T::one() - n);

        let mut result = ellippi_vc(nn, m, nm1);
        // Split calculations to avoid overflow/underflow
        result = result * -n / (T::one() - n);
        result = result * (T::one() - m) / (m - n);
        result = result + _ellipk(m) * m / (m - n);
        return Ok(result);
    }

    // Compute vc = 1-n without cancellation errors
    let vc = T::one() - n;
    Ok(ellippi_vc(n, m, vc))
}

/// Unchecked version of [ellippi].
///
/// Domain: m < 1, 0 ≤ n < 1.
pub fn _ellippi<T: Float>(n: T, m: T) -> T {
    ellippi_vc(n, m, T::one() - n)
}

#[inline]
fn ellippi_vc<T: Float>(n: T, m: T, vc: T) -> T {
    let x = T::zero();
    let y = T::one() - m;
    let z = T::one();
    let p = vc;

    _elliprf(x, y, z) + n * _elliprj(x, y, z, p) / T::from(3.0).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compare_test_data;

    fn ellippi_k(inp: &[f64]) -> f64 {
        ellippi(inp[0], inp[1] * inp[1]).unwrap()
    }

    #[test]
    fn test_ellippi() {
        compare_test_data!(
            "./tests/data/boost/ellippi2_data_f64.txt",
            ellippi_k,
            4.4e-16
        );
    }
}
