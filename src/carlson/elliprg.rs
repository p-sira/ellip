/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Boost Math, see LICENSE in this directory.
 */

use num_traits::Float;
use std::mem::swap;

use crate::{crate_util::let_mut, elliprc, elliprd, elliprf, StrErr};

// Original header from Boost Math
//  Copyright (c) 2015 John Maddock
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

/// Computes RG ([symmetric elliptic integral of the second kind](https://dlmf.nist.gov/19.16.E2_5)).
/// ```text
///                     ∞                                                             
///                 1  ⌠             t              ⎛   x       y       z   ⎞     
/// RG(x, y, z)  =  ─  ⎮ ────────────────────────── ⎜ ───── + ───── + ───── ⎟ dt
///                 4  ⎮   ________________________ ⎝ t + x   t + y   t + z ⎠     
///                    ⌡ ╲╱(t + x) (t + y) (t + z)                               
///                  0                                                             
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
/// ![Symmetric Elliptic Integral of the Second Kind](https://github.com/p-sira/ellip/blob/main/figures/elliprg_plot.svg?raw=true)
///
/// [Interactive Plot](https://github.com/p-sira/ellip/blob/main/figures/elliprg_plot.html)
///
/// # Related Functions
/// With c = csc²φ, r = 1/x², and kc² = 1 - m,
/// - [ellipe](crate::ellipe)(m) = 2 [elliprg](crate::elliprg)(0, kc², 1)
/// - [ellipeinc](crate::ellipeinc)(φ, m) = 2 [elliprg](crate::elliprg)(c - 1, c - m, c) - (c - 1) [elliprf](crate::elliprf)(c - 1, c - m, c) - [sqrt](Float::sqrt)((c - 1) * (c - m) / c)
///
/// # Examples
/// ```
/// use ellip::{elliprg, util::assert_close};
///
/// assert_close(elliprg(1.0, 0.5, 0.25).unwrap(), 0.7526721491833781, 1e-15);
/// ```
///
/// # References
/// - Maddock, John, Paul Bristow, Hubert Holin, and Xiaogang Zhang. “Boost Math Library: Special Functions - Elliptic Integrals.” Accessed April 17, 2025. <https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/ellint.html>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn elliprg<T: Float>(x: T, y: T, z: T) -> Result<T, StrErr> {
    if x < 0.0 || y < 0.0 || z < 0.0 {
        return Err("elliprg: x, y, and z must be non-negative.");
    }

    let_mut!(x, y, z);
    if x < y {
        swap(&mut x, &mut y);
    }
    if x < z {
        swap(&mut x, &mut z);
    }
    if y > z {
        swap(&mut y, &mut z);
    }

    if x == z {
        if y == z {
            return Ok(x.sqrt());
        }

        if y == 0.0 {
            return Ok(pi!() * x.sqrt() / 4.0);
        }

        if x == 0.0 {
            return Ok(y.sqrt() / 2.0);
        }

        return Ok((x * elliprc(y, x)? + y.sqrt()) / 2.0);
    }

    if y == z {
        if y == 0.0 {
            return Ok(x.sqrt() / 2.0);
        }

        return Ok((y * elliprc(x, y)? + x.sqrt()) / 2.0);
    }

    if y == 0.0 {
        let mut xn = x.sqrt();
        let mut yn = z.sqrt();
        let x0 = xn;
        let y0 = yn;
        let mut sum = 0.0;
        let mut sum_pow = 0.25;

        while (xn - yn).abs() >= 2.7 * epsilon!() * xn.abs() {
            let t = (xn * yn).sqrt();
            xn = (xn + yn) / 2.0;
            yn = t;
            sum_pow = sum_pow * 2.0;
            sum = sum + sum_pow * (xn - yn) * (xn - yn);
        }
        let rf = pi!() / (xn + yn);
        return Ok(((x0 + y0) * (x0 + y0) / 4.0 - sum) * rf / 2.0);
    }

    Ok(
        (z * elliprf(x, y, z)? - (x - z) * (y - z) * elliprd(x, y, z)? / 3.0 + (x * y / z).sqrt())
            / 2.0,
    )
}

#[cfg(not(feature = "reduce-iteration"))]
#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use super::*;
    use crate::{assert_close, compare_test_data_boost};

    fn __elliprg(inp: &[&f64]) -> f64 {
        elliprg(*inp[0], *inp[1], *inp[2]).unwrap()
    }

    fn _elliprg(inp: &[f64]) -> f64 {
        let res = elliprg(inp[0], inp[1], inp[2]).unwrap();
        inp.iter().permutations(inp.len()).skip(1).for_each(|perm| {
            assert_close!(res, __elliprg(&perm), 6.5e-16);
        });
        res
    }

    #[test]
    fn test_elliprg() {
        compare_test_data_boost!("elliprg_data.txt", _elliprg, 8.1e-16);
    }

    #[test]
    fn test_elliprg_xxx() {
        compare_test_data_boost!("elliprg_xxx.txt", _elliprg, 2.4e-16);
    }

    #[test]
    fn test_elliprg_xy0() {
        compare_test_data_boost!("elliprg_xy0.txt", _elliprg, 4.4e-16);
    }

    #[test]
    fn test_elliprg_xyy() {
        compare_test_data_boost!("elliprg_xyy.txt", _elliprg, 5.4e-16);
    }

    #[test]
    fn test_elliprg_00x() {
        compare_test_data_boost!("elliprg_00x.txt", _elliprg, f64::EPSILON);
    }

    #[test]
    fn test_elliprg_err() {
        // negative argument
        assert!(elliprg(-1.0, 1.0, 1.0).is_err());
        assert!(elliprg(1.0, -1.0, 1.0).is_err());
        assert!(elliprg(1.0, 1.0, -1.0).is_err());
    }
}
