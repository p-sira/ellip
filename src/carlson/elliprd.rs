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

use std::mem::swap;

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
    if x.min(y) < zero!() || x + y == zero!() {
        return Err("elliprd: x and y must be non-negative, and at most one can be zero.");
    }
    if z <= zero!() {
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
            return Ok(one!() / (x * x.sqrt()));
        }
        if x == zero!() {
            return Ok(three!() * pi!() / (four!() * y * y.sqrt()));
        }
    }

    if x == zero!() {
        let x0 = y.sqrt();
        let y0 = z.sqrt();
        let mut xn = x0;
        let mut yn = y0;
        let mut sum = zero!();
        let mut sum_pow = num!(0.25);

        while (xn - yn).abs() >= num!(2.7) * epsilon!() * xn.abs() {
            let t = (xn * yn).sqrt();
            xn = (xn + yn) / two!();
            yn = t;
            sum_pow = sum_pow * two!();
            let temp = xn - yn;
            sum = sum + sum_pow * temp * temp;
        }
        let rf = pi!() / (xn + yn);
        let pt = (x0 + three!() * y0) / (four!() * z * (x0 + y0)) - sum / (z * (y - z));
        return Ok(pt * rf * three!());
    }

    let mut xn = x;
    let mut yn = y;
    let mut zn = z;

    let mut an = (x + y + three!() * z) / five!();
    let a0 = an;
    let mut q = (epsilon!() / four!()).powf(-one!() / eight!())
        * (an - x).max(an - y).max(an - z)
        * num!(1.2);

    let mut fn_val = one!();
    let mut rd_sum = zero!();

    for _ in 0..N_MAX_ITERATIONS {
        let rx = xn.sqrt();
        let ry = yn.sqrt();
        let rz = zn.sqrt();
        let lambda = rx * ry + rx * rz + ry * rz;
        rd_sum = rd_sum + fn_val / (rz * (zn + lambda));
        an = (an + lambda) / four!();
        xn = (xn + lambda) / four!();
        yn = (yn + lambda) / four!();
        zn = (zn + lambda) / four!();
        fn_val = fn_val / four!();
        q = q / four!();
        if q < an {
            let x = fn_val * (a0 - x) / an;
            let y = fn_val * (a0 - y) / an;
            let z = -(x + y) / three!();
            let xyz = x * y * z;
            let z2 = z * z;
            let z3 = z2 * z;

            let e2 = x * y - six!() * z2;
            let e3 = three!() * xyz - eight!() * z3;
            let e4 = three!() * (xyz - z3) * z;
            let e5 = xyz * z2;

            let result = fn_val
                * an.powf(-num!(1.5))
                * (one!() - three!() * e2 / num!(14.0)
                    + e3 / six!()
                    + nine!() * e2 * e2 / num!(88.0)
                    - three!() * e4 / num!(22.0)
                    - nine!() * e2 * e3 / num!(52.0)
                    + three!() * e5 / num!(26.0)
                    - e2 * e2 * e2 / num!(16.0)
                    + three!() * e3 * e3 / num!(40.0)
                    + three!() * e2 * e4 / num!(20.0)
                    + num!(45.0) * e2 * e2 * e3 / num!(272.0)
                    - nine!() * (e3 * e4 + e2 * e5) / num!(68.0))
                + three!() * rd_sum;

            return Ok(result);
        }
    }

    Err("elliprd: Failed to converge.")
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
