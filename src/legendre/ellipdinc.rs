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

use crate::{ellipd, elliprd};
use num_traits::Float;

/// Compute [incomplete elliptic integral of Legendre's type](https://dlmf.nist.gov/19.2.E6).
/// ```text
///              φ
///             ⌠       sin²θ dθ
/// D(φ, m)  =  │  _________________
///             │     _____________
///             ⌡   \╱ 1 - m sin²θ
///            0
/// where m < 1
/// ```
///
/// Note that some mathematical references use the parameter k for the function,
/// where k² = m.
/// 
/// # Examples
/// ```
/// use ellip::{ellipdinc, util::assert_close};
/// use std::f64::consts::FRAC_PI_4;
/// 
/// assert_close(ellipdinc(FRAC_PI_4, 0.5).unwrap(), 0.15566274414316758, 1e-15);
/// ```
pub fn ellipdinc<T: Float>(phi: T, m: T) -> Result<T, &'static str> {
    let sign = if phi < zero!() { -one!() } else { one!() };

    let phi = phi.abs();

    if phi >= max_val!() {
        // Need to handle infinity as a special case:
        return Ok(inf!());
    }

    if phi > one!() / epsilon!() {
        // Phi is so large that phi%pi is necessarily zero (or garbage),
        // just return the second part of the duplication formula:
        return Ok(sign * two!() * phi * ellipd(m)? / pi!());
    }

    // Carlson's algorithm works only for |phi| <= pi/2,
    // use the integrand's periodicity to normalize phi
    //
    let mut rphi = phi % pi_2!();
    let mut mm = ((phi - rphi) / pi_2!()).round();
    let mut s = one!();
    if mm % two!() > half!() {
        mm = mm + one!();
        s = -one!();
        rphi = pi_2!() - rphi;
    }

    let sinp = rphi.sin();
    let cosp = rphi.cos();
    let sinp2 = sinp * sinp;
    let cosp2 = cosp * cosp;

    let c = one!() / sinp2;
    let cm1 = cosp2 / sinp2; // c - 1

    if m * sinp2 > one!() {
        return Err("ellipdinc: The argument must satisfy m sin²φ < 1.");
    }

    let mut result = zero!();
    if rphi != zero!() {
        result = s * elliprd(cm1, c - m, c)? / three!();
    }
    if mm != zero!() {
        result = result + mm * ellipd(m)?;
    }

    Ok(sign * result)
}

#[cfg(test)]
mod tests {
    use core::f64;

    use super::*;
    use crate::compare_test_data;

    fn ellipdinc_k(inp: &[f64]) -> f64 {
        ellipdinc(inp[0], inp[1] * inp[1]).unwrap()
    }

    #[test]
    fn test_ellipd() {
        compare_test_data!(
            "./tests/data/boost/ellipdinc_data.txt",
            ellipdinc_k,
            6.4e-16
        );
    }
}
