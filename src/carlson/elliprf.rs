/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Boost Math.
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
//
// Boost Software License - Version 1.0 - August 17th, 2003
//
// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:
//
// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated by
// a source language processor.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

use std::{
    f64::consts::{FRAC_PI_2, PI},
    mem::swap,
};

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
pub fn elliprf<T: Float>(x: T, y: T, z: T) -> Result<T, &'static str> {
    if x.min(y).min(z) < T::zero() || (y + z).min(x + y).min(x + z) < T::zero() {
        return Err("elliprf: x, y, and z must be non-negative, and at most one can be zero.");
    }

    // Special cases from http://dlmf.nist.gov/19.20#i
    if x == y {
        if x == z {
            // RF(x,x,x)
            return Ok(T::one() / x.sqrt());
        }

        if z == T::zero() {
            // RF(x,x,0)
            // RF(0,y,y)
            return Ok(T::from(FRAC_PI_2).unwrap() * x.sqrt());
        }

        // RF(x,x,z)
        // RF(x,y,y)
        return elliprc(z, x);
    }

    if x == z {
        if y == T::zero() {
            // RF(x,0,x)
            // RF(0,y,y)
            return Ok(T::from(PI).unwrap() / (T::from(2.0).unwrap() * x.sqrt()));
        }

        // RF(x,y,x)
        // RF(x,y,y)
        return elliprc(y, x);
    }

    if y == z {
        if x == T::zero() {
            // RF(0,y,y)
            return Ok(T::from(PI).unwrap() / (T::from(2.0).unwrap() * y.sqrt()));
        }

        // RF(x,y,y)
        return elliprc(x, y);
    }

    let mut xn = x;
    let mut yn = y;
    let mut zn = z;

    if xn == T::zero() {
        swap(&mut xn, &mut zn);
    } else if yn == T::zero() {
        swap(&mut yn, &mut zn);
    }

    if zn == T::zero() {
        let mut xn = xn.sqrt();
        let mut yn = yn.sqrt();

        while (xn - yn).abs() >= T::from(2.7).unwrap() * T::epsilon() * xn.abs() {
            let t = (xn * yn).sqrt();
            xn = (xn + yn) / T::from(2.0).unwrap();
            yn = t;
        }
        return Ok(T::from(PI).unwrap() / (xn + yn));
    }

    let res = _elliprf(xn, yn, zn);
    if res.is_nan() {
        return Err("elliprf: Failed to converge.");
    }
    Ok(res)
}

/// Unchecked version of [elliprf].
///
/// Return NAN when it fails to converge.
pub fn _elliprf<T: Float>(x: T, y: T, z: T) -> T {
    let four = T::from(4.0).unwrap();

    let mut xn = x;
    let mut yn = y;
    let mut zn = z;

    let mut an = (xn + yn + zn) / T::from(3.0).unwrap();
    let a0 = an;
    let mut q = (T::from(3.0).unwrap() * T::epsilon()).powf(-T::one() / T::from(8.0).unwrap())
        * an.abs()
            .max((an - xn).abs())
            .max((an - yn).abs())
            .max((an - zn).abs());
    let mut fn_val = T::one();
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

            return (T::one()
                + e3 * (T::from(1.0 / 14.0).unwrap()
                    + T::from(3.0).unwrap() * e3 / T::from(104.0).unwrap())
                + e2 * (T::from(-0.1).unwrap() + e2 / T::from(24.0).unwrap()
                    - (T::from(3.0).unwrap() * e3) / T::from(44.0).unwrap()
                    - T::from(5.0).unwrap() * e2 * e2 / T::from(208.0).unwrap()
                    + e2 * e3 / T::from(16.0).unwrap()))
                / an.sqrt();
        }
    }
    return T::nan();
}

const N_MAX_ITERATIONS: usize = 11;

#[cfg(test)]
mod test {
    use itertools::Itertools;

    use super::*;
    use crate::{assert_close, compare_test_data};

    fn __elliprf(inp: &Vec<&f64>) -> f64 {
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
        compare_test_data!("./tests/data/boost/ellint_rf_data.txt", _elliprf, 4.5e-16);
    }

    #[test]
    fn test_elliprf_xxx() {
        compare_test_data!("./tests/data/boost/ellint_rf_xxx.txt", _elliprf, 2.3e-16);
    }

    #[test]
    fn test_elliprf_xy0() {
        compare_test_data!("./tests/data/boost/ellint_rf_xy0.txt", _elliprf, 4.2e-16);
    }

    #[test]
    fn test_elliprf_xyy() {
        compare_test_data!("./tests/data/boost/ellint_rf_xyy.txt", _elliprf, 5.3e-16);
    }
}
