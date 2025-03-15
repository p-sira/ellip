/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Boost Math, see LICENSE in this directory.
 */

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

use num_traits::Float;

use crate::elliprf;

use super::ellipk::ellipk_precise;

/// Compute [incomplete elliptic integral of the first kind](https://dlmf.nist.gov/19.2.E4).
///
/// ```text
///              φ
///             ⌠          dθ
/// F(φ, m)  =  │  _________________
///             │     _____________
///             ⌡   \╱ 1 - m sin²θ
///            0
/// where 0 ≤ m sin²φ ≤ 1
/// ```
///
/// Note that some mathematical references use the parameter k for the function,
/// where k² = m.
///
/// # Examples
/// ```
/// use ellip::{ellipf, util::assert_close};
/// use std::f64::consts::FRAC_PI_4;
///
/// assert_close(ellipf(FRAC_PI_4, 0.5).unwrap(), 0.826017876249245, 1e-15);
/// ```
pub fn ellipf<T: Float>(phi: T, m: T) -> Result<T, &'static str> {
    let invert = if phi < zero!() { -one!() } else { one!() };
    let phi = phi.abs();

    if phi >= max_val!() {
        return Ok(invert * inf!());
    }

    if phi > one!() / epsilon!() {
        // Phi is so large that phi%pi is necessarily zero (or garbage),
        // just return the second part of the duplication formula:
        return Ok(invert * two!() * phi * ellipk_precise(m)? / pi!());
    }

    // Carlson's algorithm works only for |phi| <= pi/2,
    // use the integrand's periodicity to normalize phi
    let mut rphi = phi % pi_2!();
    let mut mm = ((phi - rphi) / pi_2!()).round();

    let mut s = one!();
    if mm % two!() > half!() {
        mm = mm + one!();
        s = -one!();
        rphi = pi_2!() - rphi;
    }

    let s2p = rphi.sin() * rphi.sin();
    let ms2p = m * s2p;
    if ms2p >= one!() {
        return Err("ellipf: m sin²φ must satisfy: 0 ≤ m sin²φ < 1.");
    }

    let c2p = rphi.cos() * rphi.cos();
    let mut result;
    if s2p > min_val!() {
        // Use http://dlmf.nist.gov/19.25#E5, note that
        // c-1 simplifies to cot^2(rphi) which avoids cancellation.
        // Likewise c - k^2 is the same as (c - 1) + (1 - k^2).
        //
        let c = one!() / s2p;
        let c_minus_one = c2p / s2p;
        let arg2 = if m != zero!() {
            let cross = (c / m).abs();
            if cross > num!(0.9) && cross < num!(1.1) {
                c_minus_one + one!() - m
            } else {
                c - m
            }
        } else {
            c
        };
        result = s * elliprf(c_minus_one, arg2, c)?;
    } else {
        result = s * rphi.sin();
    }
    if mm != zero!() {
        result = result + mm * ellipk_precise(m)?;
    }

    Ok(invert * result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compare_test_data;

    fn ellipf_k(inp: &[f64]) -> f64 {
        ellipf(inp[0], inp[1] * inp[1]).unwrap()
    }

    #[test]
    fn test_ellipf() {
        compare_test_data!("./tests/data/boost/ellipf_data.txt", ellipf_k, 5.1e-16);
    }
}
