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

use num_traits::Float;

use crate::elliprd;

/// Compute [complete elliptic integral of Legendre's type](https://dlmf.nist.gov/19.2.E8).
/// ```text
///           π/2
///          ⌠       sin²θ dθ
/// D(m)  =  │  _________________
///          │     _____________
///          ⌡   \╱ 1 - m sin²θ
///         0
/// where m < 1
/// ```
///
/// Note that some mathematical references use the parameter k for the function,
/// where k² = m.
///
/// # Examples
/// ```
/// use ellip::{ellipd, util::assert_close};
///
/// assert_close(ellipd(0.5).unwrap(), 1.0068615925073927, 1e-15);
/// ```
///
/// # Related Functions
/// - [ellipd](crate::ellipd)(m) = ([ellipk](crate::ellipk)(m) - [ellipe](crate::ellipe)(m)) / m
/// - [ellipd](crate::ellipd)(m) = [elliprd](crate::elliprd)(0, 1 - m, 1) / 3
///
/// # References
/// - Maddock, John, Paul Bristow, Hubert Holin, and Xiaogang Zhang. “Boost Math Library: Special Functions - Elliptic Integrals.” Accessed April 17, 2025. https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/ellint.html.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. https://dlmf.nist.gov/19.
///
pub fn ellipd<T: Float>(m: T) -> Result<T, &'static str> {
    if m > one!() {
        return Err("ellipd: m must be less than 1.");
    }

    // https://dlmf.nist.gov/19.2.E8
    // Using this relation, D evaluates to inf at m=1.
    if m == one!() {
        return Ok(inf!());
    }

    if m.abs() <= epsilon!() {
        return Ok(pi!() / four!());
    }

    Ok(elliprd(zero!(), one!() - m, one!())? / three!())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compare_test_data;

    fn ellipd_k(inp: &[f64]) -> f64 {
        ellipd(inp[0] * inp[0]).unwrap()
    }

    #[test]
    fn test_ellipd() {
        compare_test_data!("./tests/data/boost/ellipd_data.txt", ellipd_k, 2.9e-16);
    }
}
