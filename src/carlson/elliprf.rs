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
//  Boost.Math conceptual framework better, and to handle
//  types longer than 80-bit reals.
//  Updated 2015 to use Carlson's latest methods.

use std::mem::swap;

use num_traits::Float;

use crate::{elliprc, StrErr};

/// Computes RF ([symmetric elliptic integral of the first kind](https://dlmf.nist.gov/19.16.E1)).
/// ```text
///                     ∞                              
///                 1  ⌠             dt              
/// RF(x, y, z)  =  ─  ⎮ ───────────────────────────
///                 2  ⎮   ________________________
///                    ⌡ ╲╱(t + x) (t + y) (t + z)
///                   0                              
/// ```
///
/// ## Parameters
/// - x ∈ ℝ, x ≥ 0
/// - y ∈ ℝ, y ≥ 0
/// - z ∈ ℝ, z ≥ 0
///
/// The parameters x, y, and z are symmetric. This means swapping them does not change the value of the function.
/// At most one of them can be zero.
///
/// ## Domain
/// - Returns error if any of x, y, or z is negative, or more than one of them are zero.
///
/// ## Graph
/// ![Symmetric Elliptic Integral of the First Kind](https://github.com/p-sira/ellip/blob/main/figures/elliprf_plot.svg?raw=true)
///
/// [Interactive Plot](https://github.com/p-sira/ellip/blob/main/figures/elliprf_plot.html)
///
/// # Related Functions
/// With c = csc²φ, r = 1/x², and kc² = 1 - m,
/// - [ellipf](crate::ellipf)(φ,m) = [elliprf](crate::elliprf)(c - 1, c - m, c)
/// - [el1](crate::el1)(x, kc) = [elliprf](crate::elliprf)(r, r + m, r + 1)
///
/// # Examples
/// ```
/// use ellip::{elliprf, util::assert_close};
///
/// assert_close(elliprf(1.0, 0.5, 0.25).unwrap(), 1.370171633266872, 1e-15);
/// ```
///
/// # References
/// - Maddock, John, Paul Bristow, Hubert Holin, and Xiaogang Zhang. “Boost Math Library: Special Functions - Elliptic Integrals.” Accessed April 17, 2025. <https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/ellint.html>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn elliprf<T: Float>(x: T, y: T, z: T) -> Result<T, StrErr> {
    if x.min(y).min(z) < 0.0 || (y + z).min(x + y).min(x + z) < 0.0 {
        return Err("elliprf: x, y, and z must be non-negative, and at most one can be zero.");
    }

    // Special cases from http://dlmf.nist.gov/19.20#i
    if x == y {
        if x == z {
            // RF(x,x,x)
            return Ok(1.0 / x.sqrt());
        }

        if z == 0.0 {
            // RF(x,x,0)
            // RF(0,y,y)
            return Ok(pi!() / (2.0 * x.sqrt()));
        }

        // RF(x,x,z)
        // RF(x,y,y)
        return elliprc(z, x);
    }

    if x == z {
        if y == 0.0 {
            // RF(x,0,x)
            // RF(0,y,y)
            return Ok(pi!() / (2.0 * x.sqrt()));
        }

        // RF(x,y,x)
        // RF(x,y,y)
        return elliprc(y, x);
    }

    if y == z {
        if x == 0.0 {
            // RF(0,y,y)
            return Ok(pi!() / (2.0 * y.sqrt()));
        }

        // RF(x,y,y)
        return elliprc(x, y);
    }

    let mut xn = x;
    let mut yn = y;
    let mut zn = z;

    if xn == 0.0 {
        swap(&mut xn, &mut zn);
    } else if yn == 0.0 {
        swap(&mut yn, &mut zn);
    }

    if zn == 0.0 {
        let mut xn = xn.sqrt();
        let mut yn = yn.sqrt();

        while (xn - yn).abs() >= 2.7 * epsilon!() * xn.abs() {
            let t = (xn * yn).sqrt();
            xn = (xn + yn) / 2.0;
            yn = t;
        }
        return Ok(pi!() / (xn + yn));
    }

    let four = 4.0;

    let mut an = (xn + yn + zn) / 3.0;
    let a0 = an;
    let mut q = (3.0 * epsilon!()).powf(-1.0 / 8.0)
        * an.abs()
            .max((an - xn).abs())
            .max((an - yn).abs())
            .max((an - zn).abs());
    let mut fn_val = 1.0;
    for _ in 0..N_MAX_ITERATIONS {
        let root_x = xn.sqrt();
        let root_y = yn.sqrt();
        let root_z = zn.sqrt();

        let lambda = root_x * root_y + root_x * root_z + root_y * root_z;

        an = (an + lambda) / four;
        xn = (xn + lambda) / four;
        yn = (yn + lambda) / four;
        zn = (zn + lambda) / four;

        q = q / four;
        fn_val = fn_val * four;

        if q < an.abs() {
            let x = (a0 - x) / (an * fn_val);
            let y = (a0 - y) / (an * fn_val);
            let z = -x - y;

            let e2 = x * y - z * z;
            let e3 = x * y * z;

            return Ok((1.0
                + e3 * (1.0 / 14.0 + 3.0 * e3 / 104.0)
                + e2 * (-0.1 + e2 / 24.0 - (3.0 * e3) / 44.0 - 5.0 * e2 * e2 / 208.0
                    + e2 * e3 / 16.0))
                / an.sqrt());
        }
    }

    Err("elliprf: Failed to converge.")
}

#[cfg(not(feature = "reduce-iteration"))]
const N_MAX_ITERATIONS: usize = 11;

#[cfg(feature = "reduce-iteration")]
const N_MAX_ITERATIONS: usize = 1;

#[cfg(not(feature = "reduce-iteration"))]
#[cfg(test)]
mod tests {
    use core::f64;

    use itertools::Itertools;

    use super::*;
    use crate::{assert_close, compare_test_data_boost};

    fn __elliprf(inp: &[&f64]) -> f64 {
        elliprf(*inp[0], *inp[1], *inp[2]).unwrap()
    }

    fn _elliprf(inp: &[f64]) -> f64 {
        let res = elliprf(inp[0], inp[1], inp[2]).unwrap();
        inp.iter().permutations(inp.len()).skip(1).for_each(|perm| {
            assert_close!(res, __elliprf(&perm), 6.5e-16);
        });
        res
    }

    #[test]
    fn test_elliprf() {
        compare_test_data_boost!("elliprf_data.txt", _elliprf, 4.5e-16);
    }

    #[test]
    fn test_elliprf_xxx() {
        compare_test_data_boost!("elliprf_xxx.txt", _elliprf, 2.3e-16);
    }

    #[test]
    fn test_elliprf_xy0() {
        compare_test_data_boost!("elliprf_xy0.txt", _elliprf, 4.2e-16);
    }

    #[test]
    fn test_elliprf_xyy() {
        compare_test_data_boost!("elliprf_xyy.txt", _elliprf, 5.3e-16);
    }

    #[test]
    fn test_elliprf_0yy() {
        compare_test_data_boost!("elliprf_0yy.txt", _elliprf, f64::EPSILON);
    }

    #[test]
    fn test_elliprf_err() {
        // negative argument
        assert!(elliprf(-1.0, 1.0, 1.0).is_err());
        assert!(elliprf(1.0, -1.0, 1.0).is_err());
        assert!(elliprf(1.0, 1.0, -1.0).is_err());
        // more than one zero
        assert!(elliprf(0.0, 0.0, 1.0).is_err());
    }
}

#[cfg(feature = "reduce-iteration")]
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn force_fail_to_converge() {
        assert!(elliprf(0.2, 0.5, 1e300).is_err());
    }
}
