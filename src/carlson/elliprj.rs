/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Boost Math, see LICENSE in this directory.
 */

// Original header from Boost Math
//  Copyright (c) 2006 Xiaogang Zhang, 2015 John Maddock
//  Copyright (c) 2024 Matt Borland
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0.

use std::mem::swap;

use crate::{
    carlson::{elliprc_unchecked, elliprd_unchecked, elliprf_unchecked},
    crate_util::{case, check, declare, let_mut},
    StrErr,
};
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
/// ## Special Cases
/// - RJ(x, x, x, x) = 1/(x sqrt(x))
/// - RJ(x, y, z, z) = RD(x, y, z)
/// - RJ(x, x, x, p) = 3/(x-p) * (RC(x, p) - 1/sqrt(x)) for x ≠ p and xp ≠ 0
/// - RJ(x, y, y, y) = RD(x, y, y)
/// - RJ(x, y, y, p) = 3/(p-y) * (RC(x, y) - RC(x, p)) for y ≠ p
/// - RJ(x, y, z, p) = 0 for x = ∞ or y = ∞ or z = ∞ or p = ∞
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
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn elliprj<T: Float>(x: T, y: T, z: T, p: T) -> Result<T, StrErr> {
    check!(@neg, elliprj, "x, y, and z must be non-negative.", [x, y, z]);
    check!(@multi_zero, elliprj, [x, y, z]);

    let ans = elliprj_unchecked(x, y, z, p);

    if ans.is_finite() {
        return Ok(ans);
    }
    check!(@nan, elliprj, [x, y, z, p]);
    check!(@zero, elliprj, [p]);
    case!(@any [x, y, z, p] == inf!(), T::zero());
    Err("elliprj: Failed to converge.")
}

/// Calculate RC(1, 1 + x)
#[inline]
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
fn elliprc1p<T: Float>(y: T) -> T {
    // We can skip y = -1 check since the call from elliprj already did the check.
    // for 1 + y < 0, the integral is singular, return Cauchy principal value
    if y > 0.0 {
        y.sqrt().atan() / y.sqrt()
    } else if y == 0.0 {
        1.0
    } else if y > -0.5 {
        let arg = (-y).sqrt();
        (arg.ln_1p() - (-arg).ln_1p()) / (2.0 * (-y).sqrt())
    } else if y < -1.0 {
        (1.0 / -y).sqrt() * elliprc_unchecked(-y, -1.0 - y)
    } else {
        ((1.0 + (-y).sqrt()) / (1.0 + y).sqrt()).ln() / (-y).sqrt()
    }
}

#[inline]
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn elliprj_unchecked<T: Float>(x: T, y: T, z: T, p: T) -> T {
    let_mut!(x, y, z);
    // for p < 0, the integral is singular, return Cauchy principal value
    if p <= 0.0 {
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
        let mut value = (p - z) * elliprj_unchecked(x, y, z, p);
        value = value - 3.0 * elliprf_unchecked(x, y, z);
        value = value
            + 3.0
                * ((x * y * z) / (x * y + p * q)).sqrt()
                * elliprc_unchecked(x * y + p * q, p * q);
        return value / (z + q);
    }

    // Special cases
    // https://dlmf.nist.gov/19.20#iii
    if x == y {
        if x == z {
            if x == p {
                // RJ(x,x,x,x)
                return 1.0 / (x * x.sqrt());
            } else {
                // RJ(x,x,x,p)
                return (3.0 / (x - p)) * (elliprc_unchecked(x, p) - 1.0 / x.sqrt());
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
            return elliprd_unchecked(x, y, y);
        }
        // This prevents division by zero.
        if p.max(y) / p.min(y) > 1.2 {
            // RJ(x,y,y,p)
            return (3.0 / (p - y)) * (elliprc_unchecked(x, y) - elliprc_unchecked(x, p));
        }
    }

    if z == p {
        // RJ(x,y,z,z)
        return elliprd_unchecked(x, y, z);
    }

    declare!(mut [xn = x, yn = y, zn = z, pn = p]);
    let mut an = (x + y + z + 2.0 * p) / 5.0;
    let a0 = an;
    let mut delta = (p - x) * (p - y) * (p - z);
    let q = (epsilon!() / 5.0).powf(-1.0 / 8.0)
        * (an - x)
            .abs()
            .max((an - y).abs())
            .max((an - z).abs())
            .max((an - p).abs());

    let mut fmn = 1.0;
    let mut rc_sum = 0.0;
    let mut ans = nan!();
    for _ in 0..N_MAX_ITERATION {
        let rx = xn.sqrt();
        let ry = yn.sqrt();
        let rz = zn.sqrt();
        let rp = pn.sqrt();
        let dn = (rp + rx) * (rp + ry) * (rp + rz);
        let en = delta / (dn * dn);

        if en < -0.5 && en > -1.5 {
            let b = 2.0 * rp * (pn + rx * (ry + rz) + ry * rz) / dn;
            rc_sum = rc_sum + fmn / dn * elliprc_unchecked(1.0, b);
        } else {
            rc_sum = rc_sum + fmn / dn * elliprc1p(en);
        }

        let lambda = rx * ry + rx * rz + ry * rz;
        an = (an + lambda) / 4.0;
        fmn = fmn / 4.0;
        if fmn * q < an {
            // Calculate and return
            let x = fmn * (a0 - x) / an;
            let y = fmn * (a0 - y) / an;
            let z = fmn * (a0 - z) / an;
            let p = (-x - y - z) / 2.0;
            let xyz = x * y * z;
            let p2 = p * p;
            let p3 = p2 * p;

            let e2 = x * y + x * z + y * z - 3.0 * p2;
            let e3 = xyz + 2.0 * e2 * p + 4.0 * p3;
            let e4 = (2.0 * xyz + e2 * p + 3.0 * p3) * p;
            let e5 = xyz * p2;

            let result = fmn
                * an.powf(-1.5)
                * (1.0 - 3.0 * e2 / 14.0 + e3 / 6.0 + 9.0 * e2 * e2 / 88.0
                    - 3.0 * e4 / 22.0
                    - 9.0 * e2 * e3 / 52.0
                    + 3.0 * e5 / 26.0
                    - e2 * e2 * e2 / 16.0
                    + 3.0 * e3 * e3 / 40.0
                    + 3.0 * e2 * e4 / 20.0
                    + 45.0 * e2 * e2 * e3 / 272.0
                    - 9.0 * (e3 * e4 + e2 * e5) / 68.0);

            ans = result + 6.0 * rc_sum;
            break;
        }

        xn = (xn + lambda) / 4.0;
        yn = (yn + lambda) / 4.0;
        zn = (zn + lambda) / 4.0;
        pn = (pn + lambda) / 4.0;
        delta = delta / 64.0;
    }

    ans
}

#[cfg(not(feature = "reduce-iteration"))]
const N_MAX_ITERATION: usize = 100;

#[cfg(feature = "reduce-iteration")]
const N_MAX_ITERATION: usize = 1;

#[cfg(not(feature = "reduce-iteration"))]
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
        compare_test_data_boost!("elliprj_data.txt", _elliprj, 2.7e-14, 5e-25);
    }

    #[test]
    fn test_elliprj_e2() {
        compare_test_data_boost!("elliprj_e2.txt", _elliprj, 4.8e-14, 5e-25);
    }

    #[test]
    fn test_elliprj_e3() {
        compare_test_data_boost!("elliprj_e3.txt", _elliprj, 3.1e-15, 5e-25);
    }

    #[test]
    fn test_elliprj_e4() {
        compare_test_data_boost!("elliprj_e4.txt", _elliprj, 2.2e-16, 5e-25);
    }

    #[test]
    fn test_elliprj_zp() {
        compare_test_data_boost!("elliprj_zp.txt", _elliprj, 3.5e-15, 5e-25);
    }

    #[test]
    fn test_elliprj_special_cases() {
        use std::f64::{INFINITY, NAN};
        // Negative arguments: should return Err
        assert_eq!(
            elliprj(-1.0, 1.0, 1.0, 1.0),
            Err("elliprj: x, y, and z must be non-negative.")
        );
        assert_eq!(
            elliprj(1.0, -1.0, 1.0, 1.0),
            Err("elliprj: x, y, and z must be non-negative.")
        );
        assert_eq!(
            elliprj(1.0, 1.0, -1.0, 1.0),
            Err("elliprj: x, y, and z must be non-negative.")
        );
        // More than one zero among x, y, z: should return Err
        assert_eq!(
            elliprj(0.0, 0.0, 1.0, 1.0),
            Err("elliprj: At most one argument can be zero.")
        );
        assert_eq!(
            elliprj(0.0, 1.0, 0.0, 1.0),
            Err("elliprj: At most one argument can be zero.")
        );
        assert_eq!(
            elliprj(1.0, 0.0, 0.0, 1.0),
            Err("elliprj: At most one argument can be zero.")
        );
        // p = 0: should return Err
        assert_eq!(
            elliprj(1.0, 1.0, 1.0, 0.0),
            Err("elliprj: p cannot be zero.")
        );
        // NANs: should return Err
        assert_eq!(
            elliprj(NAN, 1.0, 1.0, 1.0),
            Err("elliprj: Arguments cannot be NAN.")
        );
        assert_eq!(
            elliprj(1.0, NAN, 1.0, 1.0),
            Err("elliprj: Arguments cannot be NAN.")
        );
        assert_eq!(
            elliprj(1.0, 1.0, NAN, 1.0),
            Err("elliprj: Arguments cannot be NAN.")
        );
        assert_eq!(
            elliprj(1.0, 1.0, 1.0, NAN),
            Err("elliprj: Arguments cannot be NAN.")
        );
        // Infs: should return zero
        assert_eq!(elliprj(INFINITY, 1.0, 1.0, 1.0).unwrap(), 0.0);
        assert_eq!(elliprj(1.0, INFINITY, 1.0, 1.0).unwrap(), 0.0);
        assert_eq!(elliprj(1.0, 1.0, INFINITY, 1.0).unwrap(), 0.0);
        assert_eq!(elliprj(1.0, 1.0, 1.0, INFINITY).unwrap(), 0.0);
    }

    #[test]
    fn cover_deadcode() {
        // else branch
        assert!(elliprc1p(-0.6).is_finite());
        // y < -1
        assert!(elliprc1p(-1.1).is_finite());
    }
}

#[cfg(feature = "reduce-iteration")]
crate::test_force_unreachable! {
    assert_eq!(elliprj(0.2, 0.5, 1e300, 1.0), Err("elliprj: Failed to converge."));
}
