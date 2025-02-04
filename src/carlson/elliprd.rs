/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Boost Math, see LICENSE in this directory.
 */

// Original header from Boost Math
//  Copyright (c) 2006 Xiaogang Zhang, 2015 John Maddock.
//  Copyright (c) 2024 Matt Borland
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  History:
//  XZ wrote the original of this file as part of the Google
//  Summer of Code 2006.  JM modified it slightly to fit into the
//  Boost.Math conceptual framework better.
//  Updated 2015 to use Carlson's latest methods.

use std::{f64::consts::PI, mem::swap};

use num_traits::Float;

/// Compute [degenerate symmetric elliptic integral of the third kind](https://dlmf.nist.gov/19.16.E5).
/// ```text
///                     ∞                                        
///                 3  ⌠                   dt                   
/// RD(x, y, z)  =  ─  ⎮ ───────────────────────────────────────
///                 2  ⎮             ___________________________
///                    ⌡ (t + z) ⋅ ╲╱(t + x) ⋅ (t + y) ⋅ (t + z)
///                   0                                        
/// where x ≥ 0, y ≥ 0, and at most one can be zero. z > 0.
/// ```
///
pub fn elliprd<T: Float>(x: T, y: T, z: T) -> Result<T, &'static str> {
    if x.min(y) < T::zero() || x + y == T::zero() {
        return Err("elliprd: x and y must be non-negative, and at most one can be zero.");
    }
    if z <= T::zero() {
        return Err("elliprd: z must be positive");
    }

    let mut x = x;
    let mut y = y;

    // Special cases
    if x == z {
        swap(&mut x, &mut y);
    }
    if y == z {
        if x == y {
            return Ok(T::one() / (x * x.sqrt()));
        }
        if x == T::zero() {
            return Ok(T::from(3.0).unwrap() * T::from(PI).unwrap()
                / (T::from(4.0).unwrap() * y * y.sqrt()));
        }
    }
    if x == T::zero() {
        let x0 = y.sqrt();
        let y0 = z.sqrt();
        let mut xn = x0;
        let mut yn = y0;
        let mut sum = T::zero();
        let mut sum_pow = T::from(0.25).unwrap();

        while (xn - yn).abs() >= T::from(2.7).unwrap() * T::epsilon() * xn.abs() {
            let t = (xn * yn).sqrt();
            xn = (xn + yn) / T::from(2.0).unwrap();
            yn = t;
            sum_pow = sum_pow * T::from(2.0).unwrap();
            let temp = xn - yn;
            sum = sum + sum_pow * temp * temp;
        }
        let rf = T::from(PI).unwrap() / (xn + yn);
        let pt = (x0 + T::from(3.0).unwrap() * y0) / (T::from(4.0).unwrap() * z * (x0 + y0))
            - sum / (z * (y - z));
        return Ok(pt * rf * T::from(3.0).unwrap());
    }

    let res = _elliprd(x, y, z);
    if res.is_nan() {
        return Err("elliprd: Failed to converge.");
    }
    Ok(res)
}

/// Unchecked version of [elliprd].
///
/// Return NAN when it fails to converge.
pub fn _elliprd<T: Float>(x: T, y: T, z: T) -> T {
    let mut xn = x;
    let mut yn = y;
    let mut zn = z;
    let mut an = (x + y + T::from(3.0).unwrap() * z) / T::from(5.0).unwrap();
    let a0 = an;
    let mut q = (T::epsilon() / T::from(4.0).unwrap()).powf(-T::one() / T::from(8.0).unwrap())
        * (an - x).max(an - y).max(an - z)
        * T::from(1.2).unwrap();

    let mut fn_val = T::one();
    let mut rd_sum = T::zero();

    for _ in 0..N_MAX_ITERATIONS {
        let rx = xn.sqrt();
        let ry = yn.sqrt();
        let rz = zn.sqrt();
        let lambda = rx * ry + rx * rz + ry * rz;
        rd_sum = rd_sum + fn_val / (rz * (zn + lambda));
        an = (an + lambda) / T::from(4.0).unwrap();
        xn = (xn + lambda) / T::from(4.0).unwrap();
        yn = (yn + lambda) / T::from(4.0).unwrap();
        zn = (zn + lambda) / T::from(4.0).unwrap();
        fn_val = fn_val / T::from(4.0).unwrap();
        q = q / T::from(4.0).unwrap();
        if q < an {
            let x = fn_val * (a0 - x) / an;
            let y = fn_val * (a0 - y) / an;
            let z = -(x + y) / T::from(3.0).unwrap();
            let xyz = x * y * z;
            let z2 = z * z;
            let z3 = z2 * z;

            let e2 = x * y - T::from(6.0).unwrap() * z2;
            let e3 = T::from(3.0).unwrap() * xyz - T::from(8.0).unwrap() * z3;
            let e4 = T::from(3.0).unwrap() * (xyz - z3) * z;
            let e5 = xyz * z2;

            let result = fn_val
                * an.powf(-T::from(1.5).unwrap())
                * (T::one() - T::from(3.0).unwrap() * e2 / T::from(14.0).unwrap()
                    + e3 / T::from(6.0).unwrap()
                    + T::from(9.0).unwrap() * e2 * e2 / T::from(88.0).unwrap()
                    - T::from(3.0).unwrap() * e4 / T::from(22.0).unwrap()
                    - T::from(9.0).unwrap() * e2 * e3 / T::from(52.0).unwrap()
                    + T::from(3.0).unwrap() * e5 / T::from(26.0).unwrap()
                    - e2 * e2 * e2 / T::from(16.0).unwrap()
                    + T::from(3.0).unwrap() * e3 * e3 / T::from(40.0).unwrap()
                    + T::from(3.0).unwrap() * e2 * e4 / T::from(20.0).unwrap()
                    + T::from(45.0).unwrap() * e2 * e2 * e3 / T::from(272.0).unwrap()
                    - T::from(9.0).unwrap() * (e3 * e4 + e2 * e5) / T::from(68.0).unwrap())
                + T::from(3.0).unwrap() * rd_sum;

            return result;
        }
    }
    T::nan()
}

const N_MAX_ITERATIONS: usize = 50;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compare_test_data;

    fn _elliprd(inp: &[f64]) -> f64 {
        elliprd(inp[0], inp[1], inp[2]).unwrap()
    }

    #[test]
    fn test_elliprd() {
        compare_test_data!("./tests/data/boost/elliprd_data.txt", _elliprd, 4.8e-16);
    }

    #[test]
    fn test_elliprd_0xy() {
        compare_test_data!("./tests/data/boost/elliprd_0xy.txt", _elliprd, 5.9e-16);
    }

    #[test]
    fn test_elliprd_0yy() {
        compare_test_data!("./tests/data/boost/elliprd_0yy.txt", _elliprd, 2.6e-16);
    }
    #[test]
    fn test_elliprd_xxx() {
        compare_test_data!("./tests/data/boost/elliprd_xxx.txt", _elliprd, 2.3e-16);
    }

    #[test]
    fn test_elliprd_xxz() {
        compare_test_data!("./tests/data/boost/elliprd_xxz.txt", _elliprd, 7.9e-16);
    }

    #[test]
    fn test_elliprd_xyy() {
        compare_test_data!("./tests/data/boost/elliprd_xyy.txt", _elliprd, 3.7e-15);
    }
}
