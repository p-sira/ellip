/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Boost Math, see LICENSE in this directory.
 */

// Original header from Boost Math
//  Copyright (c) 2006 Xiaogang Zhang, 2015 John Maddock
//  Copyright (c) 2024 Matt Borland
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  History:
//  XZ wrote the original of this file as part of the Google
//  Summer of Code 2006.  JM modified it to fit into the
//  Boost.Math conceptual framework better, and to correctly
//  handle the y < 0 case.
//  Updated 2015 to use Carlson's latest methods.

use std::f64::consts::FRAC_PI_2;

use num_traits::Float;

/// Compute [degenerate symmetric elliptic integral of RF](https://dlmf.nist.gov/19.16.E6).
/// ```text
///                  ∞                  
///              1  ⌠        dt        
/// RC(x, y)  =  ─  ⎮ ─────────────────
///              2  ⎮             _____
///                 ⌡ (t + y) ⋅ ╲╱t + x
///                0                  
/// where x ≥ 0, y ≠ 0
/// ```
///
pub fn elliprc<T: Float>(x: T, y: T) -> Result<T, &'static str> {
    if x < T::zero() {
        return Err("elliprc: x must be non-negative.");
    }

    if y == T::zero() {
        return Err("elliprc: y must be non-zero.");
    }

    Ok(_elliprc(x, y))
}

/// Unchecked version of [elliprc].
///
/// Domain: x ≥ 0, y ≠ 0
pub fn _elliprc<T: Float>(x: T, y: T) -> T {
    let mut x = x;
    let mut y = y;
    let mut prefix = T::one();
    // for y < 0, the integral is singular, return Cauchy principal value
    if y < T::zero() {
        prefix = (x / (x - y)).sqrt();
        x = x - y;
        y = -y;
    }

    if x == T::zero() {
        return prefix * T::from(FRAC_PI_2).unwrap() / y.sqrt();
    }

    if x == y {
        return prefix / x.sqrt();
    }

    if y > x {
        return prefix * ((y - x) / x).sqrt().atan() / (y - x).sqrt();
    }

    if y / x > T::from(0.5).unwrap() {
        let arg = ((x - y) / x).sqrt();
        prefix * ((arg).ln_1p() - (-arg).ln_1p()) / (T::from(2.0).unwrap() * (x - y).sqrt())
    } else {
        prefix * ((x.sqrt() + (x - y).sqrt()) / y.sqrt()).ln() / (x - y).sqrt()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::compare_test_data;

    fn _elliprc(inp: &[f64]) -> f64 {
        elliprc(inp[0], inp[1]).unwrap()
    }

    #[test]
    fn test_elliprc() {
        compare_test_data!("./tests/data/boost/ellint_rc_data.txt", _elliprc, 2.2e-16);
    }
}
