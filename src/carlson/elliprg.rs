/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Boost Math, see LICENSE in this directory.
 */

use num_traits::Float;
use std::{f64::consts::PI, mem::swap};

use crate::elliprc;

use super::{elliprd::_elliprd, elliprf::_elliprf};

// Original header from Boost Math
//  Copyright (c) 2015 John Maddock
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

/// Compute [symmetric elliptic integral of the second kind](https://dlmf.nist.gov/19.16.E2_5).
/// ```text
///                     ∞                                                             
///                 1  ⌠             t              ⎛   x       y       z   ⎞     
/// RG(x, y, z)  =  ─  ⎮ ────────────────────────── ⎜ ───── + ───── + ───── ⎟ dt
///                 4  ⎮   ________________________ ⎝ t + x   t + y   t + z ⎠     
///                    ⌡ ╲╱(t + x) (t + y) (t + z)                               
///                  0                                                             
/// where x ≥ 0, y ≥ 0, z ≥ 0
/// ```
pub fn elliprg<T: Float>(x: T, y: T, z: T) -> Result<T, &'static str> {
    if x < T::zero() || y < T::zero() || z < T::zero() {
        return Err("elliprg: x, y, and z must be non-negative.");
    }

    let mut x = x;
    let mut y = y;
    let mut z = z;
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

        if y == T::zero() {
            return Ok(T::from(PI).unwrap() * x.sqrt() / T::from(4.0).unwrap());
        }

        if x == T::zero() {
            return Ok(y.sqrt() / T::from(2.0).unwrap());
        }

        return Ok((x * elliprc(y, x)? + y.sqrt()) / T::from(2.0).unwrap());
    }

    if y == z {
        if y == T::zero() {
            return Ok(x.sqrt() / T::from(2.0).unwrap());
        }

        return Ok((y * elliprc(x, y)? + x.sqrt()) / T::from(2.0).unwrap());
    }

    if y == T::zero() {
        let mut xn = x.sqrt();
        let mut yn = z.sqrt();
        let x0 = xn;
        let y0 = yn;
        let mut sum = T::zero();
        let mut sum_pow = T::from(0.25).unwrap();

        while (xn - yn).abs() >= T::from(2.7).unwrap() * T::epsilon() * xn.abs() {
            let t = (xn * yn).sqrt();
            xn = (xn + yn) / T::from(2.0).unwrap();
            yn = t;
            sum_pow = sum_pow * T::from(2.0).unwrap();
            sum = sum + sum_pow * (xn - yn) * (xn - yn);
        }
        let rf = T::from(PI).unwrap() / (xn + yn);
        return Ok(
            ((x0 + y0) * (x0 + y0) / T::from(4.0).unwrap() - sum) * rf / T::from(2.0).unwrap()
        );
    }

    Ok(_elliprg(x, y, z))
}

/// Unchecked version of [elliprg].
pub fn _elliprg<T: Float>(x: T, y: T, z: T) -> T {
    (z * _elliprf(x, y, z) - (x - z) * (y - z) * _elliprd(x, y, z) / T::from(3.0).unwrap()
        + (x * y / z).sqrt())
        / T::from(2.0).unwrap()
}

#[cfg(test)]
mod test {
    use itertools::Itertools;

    use super::*;
    use crate::{assert_close, compare_test_data, test_util::RTOL};

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
        compare_test_data!("./tests/data/boost/ellint_rg.txt", _elliprg, 8.1e-16);
    }

    #[test]
    fn test_elliprg_xxx() {
        compare_test_data!("./tests/data/boost/ellint_rg_xxx.txt", _elliprg, 2.4e-16);
    }

    #[test]
    fn test_elliprg_xy0() {
        compare_test_data!("./tests/data/boost/ellint_rg_xy0.txt", _elliprg, 4.4e-16);
    }

    #[test]
    fn test_elliprg_xyy() {
        compare_test_data!("./tests/data/boost/ellint_rg_xyy.txt", _elliprg, 5.4e-16);
    }

    #[test]
    fn test_elliprg_00x() {
        compare_test_data!("./tests/data/boost/ellint_rg_00x.txt", _elliprg, RTOL);
    }
}
