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

use crate::elliprc;

/// Compute [symmetric elliptic integral of the first kind](https://dlmf.nist.gov/19.16.E1).
/// ```text
///                     ∞                              
///                 1  ⌠             dt              
/// RF(x, y, z)  =  ─  ⎮ ───────────────────────────
///                 2  ⎮   ________________________
///                    ⌡ ╲╱(t + x) (t + y) (t + z)
///                   0                              
/// where x ≥ 0, y ≥ 0, z ≥ 0, and at most one can be zero.
/// ```
///
/// # Examples
/// ```
/// use ellip::{elliprf, util::assert_close};
///
/// assert_close(elliprf(1.0, 0.5, 0.25).unwrap(), 1.370171633266872, 1e-15);
/// ```
pub fn elliprf<T: Float>(x: T, y: T, z: T) -> Result<T, &'static str> {
    if x.min(y).min(z) < zero!() || (y + z).min(x + y).min(x + z) < zero!() {
        return Err("elliprf: x, y, and z must be non-negative, and at most one can be zero.");
    }

    // Special cases from http://dlmf.nist.gov/19.20#i
    if x == y {
        if x == z {
            // RF(x,x,x)
            return Ok(one!() / x.sqrt());
        }

        if z == zero!() {
            // RF(x,x,0)
            // RF(0,y,y)
            return Ok(pi!() / (two!() * x.sqrt()));
        }

        // RF(x,x,z)
        // RF(x,y,y)
        return elliprc(z, x);
    }

    if x == z {
        if y == zero!() {
            // RF(x,0,x)
            // RF(0,y,y)
            return Ok(pi!() / (two!() * x.sqrt()));
        }

        // RF(x,y,x)
        // RF(x,y,y)
        return elliprc(y, x);
    }

    if y == z {
        if x == zero!() {
            // RF(0,y,y)
            return Ok(pi!() / (two!() * y.sqrt()));
        }

        // RF(x,y,y)
        return elliprc(x, y);
    }

    let mut xn = x;
    let mut yn = y;
    let mut zn = z;

    if xn == zero!() {
        swap(&mut xn, &mut zn);
    } else if yn == zero!() {
        swap(&mut yn, &mut zn);
    }

    if zn == zero!() {
        let mut xn = xn.sqrt();
        let mut yn = yn.sqrt();

        while (xn - yn).abs() >= num!(2.7) * epsilon!() * xn.abs() {
            let t = (xn * yn).sqrt();
            xn = (xn + yn) / two!();
            yn = t;
        }
        return Ok(pi!() / (xn + yn));
    }

    let four = four!();

    let mut an = (xn + yn + zn) / three!();
    let a0 = an;
    let mut q = (three!() * epsilon!()).powf(-one!() / eight!())
        * an.abs()
            .max((an - xn).abs())
            .max((an - yn).abs())
            .max((an - zn).abs());
    let mut fn_val = one!();
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

            return Ok((one!()
                + e3 * (num!(1.0 / 14.0) + three!() * e3 / num!(104.0))
                + e2 * (num!(-0.1) + e2 / num!(24.0)
                    - (three!() * e3) / num!(44.0)
                    - five!() * e2 * e2 / num!(208.0)
                    + e2 * e3 / num!(16.0)))
                / an.sqrt());
        }
    }
    Err("elliprf: Failed to converge.")
}

const N_MAX_ITERATIONS: usize = 11;

#[cfg(test)]
mod tests {
    use core::f64;

    use itertools::Itertools;

    use super::*;
    use crate::{assert_close, compare_test_data};

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
        compare_test_data!("./tests/data/boost/elliprf_data.txt", _elliprf, 4.5e-16);
    }

    #[test]
    fn test_elliprf_xxx() {
        compare_test_data!("./tests/data/boost/elliprf_xxx.txt", _elliprf, 2.3e-16);
    }

    #[test]
    fn test_elliprf_xy0() {
        compare_test_data!("./tests/data/boost/elliprf_xy0.txt", _elliprf, 4.2e-16);
    }

    #[test]
    fn test_elliprf_xyy() {
        compare_test_data!("./tests/data/boost/elliprf_xyy.txt", _elliprf, 5.3e-16);
    }

    #[test]
    fn test_elliprf_0yy() {
        compare_test_data!("./tests/data/boost/elliprf_0yy.txt", _elliprf, f64::EPSILON);
    }
}
