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
//  handle the p < 0 case.
//  Updated 2015 to use Carlson's latest methods.

use std::mem::swap;

use crate::{elliprc, elliprd, elliprf};
use num_traits::Float;

/// Computes RJ ([symmetric elliptic integral of the third kind](https://dlmf.nist.gov/19.16.E2)).
/// ```text
///                        ∞                              
///                    3  ⌠                  dt              
/// RJ(x, y, z, p)  =  ─  ⎮ ─────────────────────────────────────
///                    2  ⎮             _______________________
///                       ⌡ (t + p) ⋅ ╲╱(t + x) (t + y) (t + z)
///                      0
/// ```
///
/// ## Parameters
/// - x ∈ ℝ, x ≥ 0
/// - y ∈ ℝ, y ≥ 0
/// - z ∈ ℝ, z ≥ 0
/// - p ∈ ℝ, p ≠ 0
///
/// The parameters x, y, and z are symmetric. This means swapping them does not change the value of the function.
/// At most one of them can be zero.
///
/// ## Domain
/// - Returns error if:
///   - any of x, y, or z is negative, or more than one of them are zero,
///   - or p = 0.
/// - Returns the Cauchy principal value if p < 0.
///
/// ## Graph
/// ![Symmetric Elliptic Integral of the Third Kind](https://github.com/p-sira/ellip/blob/main/figures/elliprj_plot.svg?raw=true)
///
/// [Interactive Plot](https://github.com/p-sira/ellip/blob/main/figures/elliprj_plot.html)
///
/// # Related Functions
/// With c = csc²φ and kc² = 1 - m,
/// - [ellippi](crate::ellippi)(φ, n, m) = n / 3 * [elliprj](crate::elliprj)(c - 1, c - m, c, c - n) + [ellipf](crate::ellipf)(φ, m)
///
/// # Examples
/// ```
/// use ellip::{elliprj, util::assert_close};
///
/// assert_close(elliprj(1.0, 0.5, 0.25, 0.125).unwrap(), 5.680557292035963, 1e-15);
/// ```
///
/// # References
/// - Maddock, John, Paul Bristow, Hubert Holin, and Xiaogang Zhang. “Boost Math Library: Special Functions - Elliptic Integrals.” Accessed April 17, 2025. <https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/ellint.html>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
///
pub fn elliprj<T: Float>(x: T, y: T, z: T, p: T) -> Result<T, &'static str> {
    if x.min(y).min(z) < zero!() || (y + z).min(x + y).min(x + z) == zero!() {
        return Err("elliprj: x, y, and z must be non-negative, and at most one can be zero.");
    }
    if p == zero!() {
        return Err("elliprj: p must be non-zero");
    }

    let mut x = x;
    let mut y = y;
    let mut z = z;

    // for p < 0, the integral is singular, return Cauchy principal value
    if p < zero!() {
        // We must ensure that x < y < z.
        // Since the integral is symmetrical in x, y and z
        // we can just permute the values:
        if x > y {
            swap(&mut x, &mut y);
        }
        if y > z {
            swap(&mut y, &mut z);
        }
        if x > y {
            swap(&mut x, &mut y);
        }

        let q = -p;
        let p = (z * (x + y + q) - x * y) / (z + q);
        let mut value = (p - z) * _elliprj(x, y, z, p)?;
        value = value - three!() * elliprf(x, y, z)?;
        value = value
            + three!() * ((x * y * z) / (x * y + p * q)).sqrt() * elliprc(x * y + p * q, p * q)?;
        return Ok(value / (z + q));
    }

    _elliprj(x, y, z, p)
}

/// Calculate RC(1, 1 + x)
#[inline]
fn elliprc1p<T: Float>(y: T) -> Result<T, &'static str> {
    // We can skip this check since the call from elliprj already did the check.
    // if y == -1.0 {
    //     return Err("elliprc1p: y cannot be -1.0.");
    // }

    // for 1 + y < 0, the integral is singular, return Cauchy principal value
    if y < -one!() {
        Ok((one!() / -y).sqrt() * elliprc(-y, -one!() - y)?)
    } else if y == zero!() {
        Ok(one!())
    } else if y > zero!() {
        Ok(y.sqrt().atan() / y.sqrt())
    } else if y > num!(-0.5) {
        let arg = (-y).sqrt();
        Ok((arg.ln_1p() - (-arg).ln_1p()) / (two!() * (-y).sqrt()))
    } else {
        Ok(((one!() + (-y).sqrt()) / (one!() + y).sqrt()).ln() / (-y).sqrt())
    }
}

#[inline]
fn _elliprj<T: Float>(x: T, y: T, z: T, p: T) -> Result<T, &'static str> {
    let mut x = x;
    let mut z = z;
    // Special cases
    // https://dlmf.nist.gov/19.20#iii
    if x == y {
        if x == z {
            if x == p {
                // RJ(x,x,x,x)
                return Ok(one!() / (x * x.sqrt()));
            } else {
                // RJ(x,x,x,p)
                return Ok((three!() / (x - p)) * (elliprc(x, p)? - one!() / x.sqrt()));
            }
        } else {
            // RJ(x,x,z,p)
            swap(&mut x, &mut z);
            // RJ(x,y,y,p)
            // Fall through to next if block.
        }
    }

    if y == z {
        if y == p {
            // RJ(x,y,y,y)
            return elliprd(x, y, y);
        }
        // This prevents division by zero.
        if p.max(y) / p.min(y) > num!(1.2) {
            // RJ(x,y,y,p)
            return Ok((three!() / (p - y)) * (elliprc(x, y)? - elliprc(x, p)?));
        }
    }

    if z == p {
        // RJ(x,y,z,z)
        return elliprd(x, y, z);
    }

    let mut xn = x;
    let mut yn = y;
    let mut zn = z;
    let mut pn = p;
    let mut an = (x + y + z + two!() * p) / five!();
    let a0 = an;
    let mut delta = (p - x) * (p - y) * (p - z);
    let q = (epsilon!() / five!()).powf(-one!() / eight!())
        * (an - x)
            .abs()
            .max((an - y).abs())
            .max((an - z).abs())
            .max((an - p).abs());

    let mut fmn = one!();
    let mut rc_sum = zero!();

    for _ in 0..N_MAX_ITERATION {
        let rx = xn.sqrt();
        let ry = yn.sqrt();
        let rz = zn.sqrt();
        let rp = pn.sqrt();
        let dn = (rp + rx) * (rp + ry) * (rp + rz);
        let en = delta / (dn * dn);

        if en < num!(-0.5) && en > num!(-1.5) {
            let b = two!() * rp * (pn + rx * (ry + rz) + ry * rz) / dn;
            rc_sum = rc_sum + fmn / dn * elliprc(one!(), b)?;
        } else {
            rc_sum = rc_sum + fmn / dn * elliprc1p(en)?;
        }

        let lambda = rx * ry + rx * rz + ry * rz;
        an = (an + lambda) / four!();
        fmn = fmn / four!();

        if fmn * q < an {
            // Calculate and return
            let x = fmn * (a0 - x) / an;
            let y = fmn * (a0 - y) / an;
            let z = fmn * (a0 - z) / an;
            let p = (-x - y - z) / two!();
            let xyz = x * y * z;
            let p2 = p * p;
            let p3 = p2 * p;

            let e2 = x * y + x * z + y * z - three!() * p2;
            let e3 = xyz + two!() * e2 * p + four!() * p3;
            let e4 = (two!() * xyz + e2 * p + three!() * p3) * p;
            let e5 = xyz * p2;

            let result = fmn
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
                    - nine!() * (e3 * e4 + e2 * e5) / num!(68.0));

            return Ok(result + six!() * rc_sum);
        }

        xn = (xn + lambda) / four!();
        yn = (yn + lambda) / four!();
        zn = (zn + lambda) / four!();
        pn = (pn + lambda) / four!();
        delta = delta / num!(64.0);
    }
    Err("elliprj: Fail to converge")
}

const N_MAX_ITERATION: usize = 100;

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use super::*;
    use crate::{assert_close, compare_test_data_boost};

    fn __elliprj(inp: &[&f64]) -> f64 {
        elliprj(*inp[0], *inp[1], *inp[2], *inp[3]).unwrap()
    }

    fn _elliprj(inp: &[f64]) -> f64 {
        let res = elliprj(inp[0], inp[1], inp[2], inp[3]).unwrap();
        let (p, sym_params) = inp.split_last().unwrap();
        sym_params
            .iter()
            .permutations(sym_params.len())
            .skip(1)
            .for_each(|mut perm| {
                perm.push(p);
                assert_close!(res, __elliprj(&perm), 3e-14);
            });
        res
    }

    #[test]
    fn test_elliprj() {
        compare_test_data_boost!(
            "./tests/data/boost/elliprj_data.txt",
            _elliprj,
            2.7e-14,
            5e-25
        );
    }

    #[test]
    fn test_elliprj_e2() {
        compare_test_data_boost!(
            "./tests/data/boost/elliprj_e2.txt",
            _elliprj,
            4.8e-14,
            5e-25
        );
    }

    #[test]
    fn test_elliprj_e3() {
        compare_test_data_boost!(
            "./tests/data/boost/elliprj_e3.txt",
            _elliprj,
            3.1e-15,
            5e-25
        );
    }

    #[test]
    fn test_elliprj_e4() {
        compare_test_data_boost!(
            "./tests/data/boost/elliprj_e4.txt",
            _elliprj,
            2.2e-16,
            5e-25
        );
    }

    #[test]
    fn test_elliprj_zp() {
        compare_test_data_boost!(
            "./tests/data/boost/elliprj_zp.txt",
            _elliprj,
            3.5e-15,
            5e-25
        );
    }
}
