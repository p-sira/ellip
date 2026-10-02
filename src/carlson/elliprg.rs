/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Boost Math, see LICENSE in this directory.
 */

use core::mem::swap;
use num_traits::Float;

use crate::{
    carlson::elliprc_unchecked,
    crate_util::{check, let_mut},
    StrErr,
};

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
///
/// ## Domain
/// - Returns error if any of x, y, or z is negative or infinite.
///
/// ## Graph
/// ![Symmetric Elliptic Integral of the Second Kind](https://github.com/p-sira/ellip/blob/main/figures/elliprg.svg?raw=true)
///
/// [Interactive Plot](https://p-sira.github.io/ellippy/_static/figures/elliprg.html)
///
/// ## Special Cases
/// - RG(x, x, x) = sqrt(x)
/// - RG(0, y, y) = π/4 * sqrt(y)
/// - RG(x, y, y) = (y * RC(x, y) + sqrt(x))/2
/// - RG(0, 0, z) = sqrt(z)/2
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
pub fn elliprg<T: Float>(x: T, y: T, z: T) -> Result<T, StrErr> {
    check!(@neg, elliprg, [x, y, z]);

    let ans = elliprg_unchecked(x, y, z);
    if ans.is_finite() {
        return Ok(ans);
    }
    check!(@nan, elliprg, [x, y, z]);
    Err("elliprg: Arguments must be finite.")
}

/// Unsafe version of [elliprg](crate::elliprg).
/// <div class="warning">⚠️ Unstable feature. May subject to changes.</div>
///
/// Undefined behavior with invalid arguments and edge cases.
/// # Known Invalid Cases
/// - x < 0, y < 0, z < 0
/// - x, y, or z are infinite.
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub fn elliprg_unchecked<T: Float>(x: T, y: T, z: T) -> T {
    // Homogeneity keeps products in duplication and special cases in range.
    // Separate rescaling factors: scale*sqrt(scale) can itself overflow.
    let scale = x.abs().max(y.abs()).max(z.abs());
    let limit = T::max_value().sqrt().sqrt().sqrt();
    if scale.is_finite() && scale > 0.0 && (scale > limit || scale < 1.0 / limit)
        // Do not turn nonzero arguments into singular zeros at extreme ratios.
        && [x, y, z].iter().all(|&v| v == 0.0 || v / scale != 0.0)
    {
        return elliprg_unchecked(x / scale, y / scale, z / scale) * scale.sqrt();
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
            return x.sqrt();
        }

        if y == 0.0 {
            return pi!() * x.sqrt() / 4.0;
        }

        return (x * elliprc_unchecked(y, x) + y.sqrt()) / 2.0;
    }

    if y == z {
        if y == 0.0 {
            return x.sqrt() / 2.0;
        }

        return (y * elliprc_unchecked(x, y) + x.sqrt()) / 2.0;
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
        return ((x0 + y0) * (x0 + y0) / 4.0 - sum) * rf / 2.0;
    }

    let (rf, rd) = elliprf_rd_unchecked(x, y, z);
    (z * rf - (x - z) * (y - z) * rd / 3.0 + (x * y / z).sqrt()) / 2.0
}

#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub(crate) fn elliprf_rd_unchecked<T: Float>(x: T, y: T, z: T) -> (T, T) {
    let mut xn = x;
    let mut yn = y;
    let mut zn = z;
    let a0_rf = (xn + yn + zn) / 3.0;
    let a0_rd = (xn + yn + 3.0 * zn) / 5.0;
    let mut an_rf = a0_rf;
    let mut an_rd = a0_rd;
    let eps = epsilon!();
    let mut q_rf = (3.0 * eps).powf(-1.0 / 8.0)
        * an_rf
            .abs()
            .max((an_rf - x).abs())
            .max((an_rf - y).abs())
            .max((an_rf - z).abs());
    let mut q_rd = (eps / 4.0).powf(-1.0 / 8.0) * (an_rd - x).max(an_rd - y).max(an_rd - z) * 1.2;

    let mut fn_rf = 1.0;
    let mut fn_rd = 1.0;
    let mut rd_sum = 0.0;
    let mut rf_res = nan!();
    let mut rd_res = nan!();
    let mut rf_done = false;
    let mut rd_done = false;

    for _ in 0..N_MAX_ITERATIONS {
        let rx = xn.sqrt();
        let ry = yn.sqrt();
        let rz = zn.sqrt();
        let lambda = rx * ry + rx * rz + ry * rz;

        if !rd_done {
            rd_sum = rd_sum + fn_rd / (rz * (zn + lambda));
            an_rd = (an_rd + lambda) / 4.0;
            fn_rd = fn_rd / 4.0;
            q_rd = q_rd / 4.0;
            if q_rd < an_rd {
                let x_c = fn_rd * (a0_rd - x) / an_rd;
                let y_c = fn_rd * (a0_rd - y) / an_rd;
                let z_c = -(x_c + y_c) / 3.0;
                let xyz = x_c * y_c * z_c;
                let z2 = z_c * z_c;
                let z3 = z2 * z_c;
                let e2 = x_c * y_c - 6.0 * z2;
                let e3 = 3.0 * xyz - 8.0 * z3;
                let e4 = 3.0 * (xyz - z3) * z_c;
                let e5 = xyz * z2;
                let inv_an_1_5 = an_rd.powf(-1.5);
                rd_res = fn_rd
                    * inv_an_1_5
                    * (1.0 - 3.0 * e2 / 14.0 + e3 / 6.0 + 9.0 * e2 * e2 / 88.0
                        - 3.0 * e4 / 22.0
                        - 9.0 * e2 * e3 / 52.0
                        + 3.0 * e5 / 26.0
                        - e2 * e2 * e2 / 16.0
                        + 3.0 * e3 * e3 / 40.0
                        + 3.0 * e2 * e4 / 20.0
                        + 45.0 * e2 * e2 * e3 / 272.0
                        - 9.0 * (e3 * e4 + e2 * e5) / 68.0)
                    + 3.0 * rd_sum;
                rd_done = true;
            }
        }

        if !rf_done {
            an_rf = (an_rf + lambda) / 4.0;
            q_rf = q_rf / 4.0;
            fn_rf = fn_rf * 4.0;
            if q_rf < an_rf.abs() {
                let x_c = (a0_rf - x) / (an_rf * fn_rf);
                let y_c = (a0_rf - y) / (an_rf * fn_rf);
                let z_c = -x_c - y_c;
                let e2 = x_c * y_c - z_c * z_c;
                let e3 = x_c * y_c * z_c;
                rf_res = (1.0
                    + e3 * (1.0 / 14.0 + 3.0 * e3 / 104.0)
                    + e2 * (-0.1 + e2 / 24.0
                        - (3.0 * e3) / 44.0
                        - 5.0 * e2 * e2 / 208.0
                        + e2 * e3 / 16.0))
                    / an_rf.sqrt();
                rf_done = true;
            }
        }

        if rf_done && rd_done {
            break;
        }

        xn = (xn + lambda) / 4.0;
        yn = (yn + lambda) / 4.0;
        zn = (zn + lambda) / 4.0;
    }

    (rf_res, rd_res)
}

#[cfg(not(feature = "test_force_fail"))]
const N_MAX_ITERATIONS: usize = 50;

#[cfg(feature = "test_force_fail")]
const N_MAX_ITERATIONS: usize = 1;

#[cfg(not(feature = "test_force_fail"))]
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
    fn test_elliprg_special_cases() {
        use std::f64::{INFINITY, NAN};
        // Negative arguments: should return Err
        assert_eq!(
            elliprg(-1.0, 1.0, 1.0),
            Err("elliprg: Arguments must be non-negative.")
        );
        assert_eq!(
            elliprg(1.0, -1.0, 1.0),
            Err("elliprg: Arguments must be non-negative.")
        );
        assert_eq!(
            elliprg(1.0, 1.0, -1.0),
            Err("elliprg: Arguments must be non-negative.")
        );
        // NANs: should return Err
        assert_eq!(
            elliprg(NAN, 1.0, 1.0),
            Err("elliprg: Arguments cannot be NAN.")
        );
        assert_eq!(
            elliprg(1.0, NAN, 1.0),
            Err("elliprg: Arguments cannot be NAN.")
        );
        assert_eq!(
            elliprg(1.0, 1.0, NAN),
            Err("elliprg: Arguments cannot be NAN.")
        );
        // Infinity arguments should return Err
        assert_eq!(
            elliprg(INFINITY, 1.0, 1.0),
            Err("elliprg: Arguments must be finite.")
        );
        assert_eq!(
            elliprg(1.0, INFINITY, 1.0),
            Err("elliprg: Arguments must be finite.")
        );
        assert_eq!(
            elliprg(1.0, 1.0, INFINITY),
            Err("elliprg: Arguments must be finite.")
        );
    }

    // Regression for audit finding A8: https://github.com/p-sira/ellip/pull/122
    #[test]
    fn test_elliprg_extreme_scales() {
        for (scale, expected) in [
            (1e200_f64, 1.4018470999908951e100),
            (1e-200_f64, 1.4018470999908951e-100),
        ] {
            let actual = elliprg(scale, 2.0 * scale, 3.0 * scale).unwrap();
            assert!(actual.is_finite());
            assert!((actual - expected).abs() <= 2e-15 * expected);
        }
        for scale in [1e-200_f64, 1e-100, 1e100, 1e200] {
            let actual = elliprg(scale, 2.0 * scale, 3.0 * scale).unwrap() / scale.sqrt();
            let expected = elliprg(1.0, 2.0, 3.0).unwrap();
            assert!((actual - expected).abs() <= 3e-15 * expected.abs());
        }
        let actual = elliprg(1e30_f32, 2e30, 3e30).unwrap() as f64;
        let expected = 1.4018470999908951e15;
        assert!((actual - expected).abs() <= 5e-7 * expected);
    }
}
