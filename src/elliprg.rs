/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Boost Math.
 */

use std::f64::consts::PI;

use crate::{elliprc, elliprd, elliprf};

// Ellipe modified how the function handles special cases by swapping
// the variables to arrange them first, then check against the special
// case conditions.

// Original header from Boost Math
//  Copyright (c) 2015 John Maddock
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
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
pub fn elliprg(x: f64, y: f64, z: f64) -> Result<f64, &'static str> {
    if x < 0.0 || y < 0.0 || z < 0.0 {
        return Err("elliprg: x, y, and z must be non-negative.");
    }

    let mut x = x;
    let mut y = y;
    let mut z = z;
    // Function is symmetric in x, y and z, but we require
    // (x - z)(y - z) >= 0 to avoid cancellation error in the result
    // which implies (for example) x >= z >= y
    if x < y {
        std::mem::swap(&mut x, &mut y);
    }
    if x < z {
        std::mem::swap(&mut x, &mut z);
    }
    if y > z {
        std::mem::swap(&mut y, &mut z);
    }

    // Special cases from http://dlmf.nist.gov/19.20#ii
    if x == z {
        if y == z {
            // RG(x,x,x)
            return Ok(x.sqrt());
        }

        if y == 0.0 {
            // RG(x,0,x)
            // RG(0,y,y)
            return Ok(PI * x.sqrt() / 4.0);
        }

        if x == 0.0 {
            // RG(0,y,0)
            // RG(0,0,z)
            return Ok(y.sqrt() / 2.0);
        }

        // RG(x,y,x)
        // 2 * RG(x,y,y)
        return Ok((x * elliprc(y, x)? + y.sqrt()) / 2.0);
    }

    if y == z {
        if y == 0.0 {
            // RG(x,0,0)
            // RG(0,0,z)
            return Ok(x.sqrt() / 2.0);
        }

        // RG(x,y,y)
        return Ok((y * elliprc(x, y)? + x.sqrt()) / 2.0);
    }

    if y == 0.0 {
        // Special case for y = 0
        let mut xn = x.sqrt();
        let mut yn = z.sqrt(); // Swap z with y
        let x0 = xn;
        let y0 = yn;
        let mut sum = 0.0;
        let mut sum_pow = 0.25;

        while (xn - yn).abs() >= 2.7 * f64::EPSILON * xn.abs() {
            let t = (xn * yn).sqrt();
            xn = (xn + yn) / 2.0;
            yn = t;
            sum_pow *= 2.0;
            sum += sum_pow * (xn - yn).powi(2);
        }
        let rf = PI / (xn + yn);
        return Ok(((x0 + y0).powi(2) / 4.0 - sum) * rf / 2.0);
    }

    Ok(
        (z * elliprf(x, y, z)? - (x - z) * (y - z) * elliprd(x, y, z)? / 3.0 + (x * y / z).sqrt())
            / 2.0,
    )
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::compare_test_data;

    fn _elliprg(inp: &[f64]) -> f64 {
        elliprg(inp[0], inp[1], inp[2]).unwrap()
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
}
