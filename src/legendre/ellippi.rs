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

use crate::{ellipk, elliprf, elliprj};

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
    if m >= one!() {
        return Err("ellippi: m must be less than 1.");
    }
    if n > one!() {
        return Err("ellippi: n must be less than 1.");
    }

    if n == zero!() {
        if m == zero!() {
            return Ok(pi_2!());
        }
        return ellipk(m);
    }

    if n < zero!() {
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
fn ellippi_vc<T: Float>(n: T, m: T, vc: T) -> Result<T, &'static str> {
    let x = zero!();
    let y = one!() - m;
    let z = one!();
    let p = vc;

    Ok(elliprf(x, y, z)? + n * elliprj(x, y, z, p)? / three!())
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
            "./tests/data/boost/ellint_pi2_data_f64.txt",
            ellippi_k,
            4.4e-16
        );
    }
}
