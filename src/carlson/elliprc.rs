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

use num_traits::Float;

use crate::{crate_util::let_mut, StrErr};

/// Computes RC ([degenerate symmetric elliptic integral of RF](https://dlmf.nist.gov/19.16.E6)).
/// ```text
///                  ∞                  
///              1  ⌠        dt        
/// RC(x, y)  =  ─  ⎮ ─────────────────
///              2  ⎮             _____
///                 ⌡ (t + y) ⋅ ╲╱t + x
///                0                  
/// ```
///
/// ## Parameters
/// - x ∈ ℝ, x ≥ 0
/// - y ∈ ℝ, y ≠ 0
///
/// ## Domain
/// - Returns error if x < 0 or y = 0.
/// - Returns the Cauchy principal value if y < 0.
///
/// ## Graph
/// ![Degenerate Symmetric Elliptic Integral of RF (RC)](https://github.com/p-sira/ellip/blob/main/figures/elliprc_plot.svg?raw=true)
///
/// [Interactive Plot](https://github.com/p-sira/ellip/blob/main/figures/elliprc_plot.html)
///
/// ## Notes
/// RC is a degenerate case of the RF. It is an elementary function rather than an elliptic integral.
///
/// # Related Functions
/// - [elliprc](crate::elliprc)(x, y) = [elliprf](crate::elliprf)(x, y, y)
///
/// # Examples
/// ```
/// use ellip::{elliprc, util::assert_close};
///
/// assert_close(elliprc(1.0, 0.5).unwrap(), 1.2464504802804608, 1e-15);
/// ```
///
/// # References
/// - Maddock, John, Paul Bristow, Hubert Holin, and Xiaogang Zhang. “Boost Math Library: Special Functions - Elliptic Integrals.” Accessed April 17, 2025. <https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/ellint.html>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
/// - The SciPy Community. “SciPy: Special Functions - Elliprc.” Accessed April 17, 2025. <https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.elliprc.html>.
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn elliprc<T: Float>(x: T, y: T) -> Result<T, StrErr> {
    if x < 0.0 {
        return Err("elliprc: x must be non-negative.");
    }

    if y == 0.0 {
        return Err("elliprc: y must be non-zero.");
    }

    Ok(_elliprc(x, y))
}

#[inline]
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
fn _elliprc<T: Float>(x: T, y: T) -> T {
    let_mut!(x, y);
    let mut prefix = 1.0;
    // for y < 0, the integral is singular, return Cauchy principal value
    if y < 0.0 {
        prefix = (x / (x - y)).sqrt();
        x = x - y;
        y = -y;
    }

    if x == 0.0 {
        return prefix * pi_2!() / y.sqrt();
    }

    if x == y {
        return prefix / x.sqrt();
    }

    if y > x {
        return prefix * ((y - x) / x).sqrt().atan() / (y - x).sqrt();
    }

    if y / x > 0.5 {
        let arg = ((x - y) / x).sqrt();
        prefix * ((arg).ln_1p() - (-arg).ln_1p()) / (2.0 * (x - y).sqrt())
    } else {
        prefix * ((x.sqrt() + (x - y).sqrt()) / y.sqrt()).ln() / (x - y).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compare_test_data_boost;

    fn _elliprc(inp: &[f64]) -> f64 {
        elliprc(inp[0], inp[1]).unwrap()
    }

    #[test]
    fn test_elliprc() {
        compare_test_data_boost!("elliprc_data.txt", _elliprc, f64::EPSILON);
    }

    #[test]
    fn test_elliprc_err() {
        // x < 0
        assert!(elliprc(-1.0, 1.0).is_err());
        // y == 0
        assert!(elliprc(1.0, 0.0).is_err());
    }
}
