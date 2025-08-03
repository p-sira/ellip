/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Boost Math, see LICENSE in this directory.
 */

// Original header from Boost Math
//  Copyright (c) 2006 Xiaogang Zhang, 2015 John Maddock.
//  Copyright (c) 2024 Matt Borland
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0.

use std::mem::swap;

use num_traits::Float;

use crate::{
    crate_util::{case, check, let_mut, return_if_valid_else},
    StrErr,
};

/// Computes RD ([degenerate symmetric elliptic integral of the third kind](https://dlmf.nist.gov/19.16.E5)).
/// ```text
///                     ∞                                        
///                 3  ⌠                   dt                   
/// RD(x, y, z)  =  ─  ⎮ ───────────────────────────────────────
///                 2  ⎮             ___________________________
///                    ⌡ (t + z) ⋅ ╲╱(t + x) ⋅ (t + y) ⋅ (t + z)
///                   0                                        
/// ```
///
/// ## Parameters
/// - x ∈ ℝ, x ≥ 0
/// - y ∈ ℝ, y ≥ 0
/// - z ∈ ℝ, z > 0
///
/// The parameters x and y (but not z!) are symmetric. This means swapping them does not change the value of the function.
/// At most one of them can be zero.
///
/// ## Domain
/// - Returns error if x < 0, y < 0, z ≤ 0 or when both x and y are zero.
///
/// ## Graph
/// ![Degenerate Symmetric Elliptic Integral of the Third Kind](https://github.com/p-sira/ellip/blob/main/figures/elliprd_plot.svg?raw=true)
///
/// [Interactive Plot](https://github.com/p-sira/ellip/blob/main/figures/elliprd_plot.html)
///
/// ## Special Cases
/// - RD(x, x, x) = 1/(x sqrt(x))
/// - RD(0, y, y) = 3/4 * π / (y sqrt(y))
/// - RD(x, y, y) = 3/(2(y-x)) * (RC(x, y) - sqrt(x)/y) for x ≠ y
/// - RD(x, x, z) = 3/(z-x) * (RC(z, x) - 1/sqrt(z)) for x ≠ z
/// - RD(x, y, z) = 0 for x = ∞ or y = ∞ or z = ∞
///
/// # Related Functions
/// With c = csc²φ,
/// - [ellipdinc](crate::ellipdinc)(φ, m) = [elliprd](crate::elliprd)(c - 1, c - m, c) / 3
///
/// # Examples
/// ```
/// use ellip::{elliprd, util::assert_close};
///
/// assert_close(elliprd(1.0, 0.5, 0.25).unwrap(), 4.022594757168912, 1e-15);
/// ```
///
/// # References
/// - Maddock, John, Paul Bristow, Hubert Holin, and Xiaogang Zhang. “Boost Math Library: Special Functions - Elliptic Integrals.” Accessed April 17, 2025. <https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/ellint.html>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn elliprd<T: Float>(x: T, y: T, z: T) -> Result<T, StrErr> {
    if x.min(y) < 0.0 || x + y == 0.0 {
        return Err("elliprd: x and y must be non-negative, and at most one can be zero.");
    }
    if z <= 0.0 {
        return Err("elliprd: z must be positive");
    }

    let_mut!(x, y);

    // Special cases
    if x == z {
        swap(&mut x, &mut y);
    }

    if y == z {
        if x == y {
            return Ok(1.0 / (x * x.sqrt()));
        }
        if x == 0.0 {
            return Ok(3.0 * pi!() / (4.0 * y * y.sqrt()));
        }
    }

    if x == 0.0 {
        let x0 = y.sqrt();
        let y0 = z.sqrt();
        let mut xn = x0;
        let mut yn = y0;
        let mut sum = 0.0;
        let mut sum_pow = 0.25;

        while (xn - yn).abs() >= 2.7 * epsilon!() * xn.abs() {
            let t = (xn * yn).sqrt();
            xn = (xn + yn) / 2.0;
            yn = t;
            sum_pow = sum_pow * 2.0;
            let temp = xn - yn;
            sum = sum + sum_pow * temp * temp;
        }
        let rf = pi!() / (xn + yn);
        let pt = (x0 + 3.0 * y0) / (4.0 * z * (x0 + y0)) - sum / (z * (y - z));
        return Ok(pt * rf * 3.0);
    }

    let mut xn = x;
    let mut yn = y;
    let mut zn = z;

    let mut an = (x + y + 3.0 * z) / 5.0;
    let a0 = an;
    let mut q = (epsilon!() / 4.0).powf(-1.0 / 8.0) * (an - x).max(an - y).max(an - z) * 1.2;

    let mut fn_val = 1.0;
    let mut rd_sum = 0.0;

    let mut ans = nan!();
    for _ in 0..N_MAX_ITERATIONS {
        let rx = xn.sqrt();
        let ry = yn.sqrt();
        let rz = zn.sqrt();
        let lambda = rx * ry + rx * rz + ry * rz;
        rd_sum = rd_sum + fn_val / (rz * (zn + lambda));
        an = (an + lambda) / 4.0;
        xn = (xn + lambda) / 4.0;
        yn = (yn + lambda) / 4.0;
        zn = (zn + lambda) / 4.0;
        fn_val = fn_val / 4.0;
        q = q / 4.0;
        if q < an {
            let x = fn_val * (a0 - x) / an;
            let y = fn_val * (a0 - y) / an;
            let z = -(x + y) / 3.0;
            let xyz = x * y * z;
            let z2 = z * z;
            let z3 = z2 * z;

            let e2 = x * y - 6.0 * z2;
            let e3 = 3.0 * xyz - 8.0 * z3;
            let e4 = 3.0 * (xyz - z3) * z;
            let e5 = xyz * z2;

            ans = fn_val
                * an.powf(-1.5)
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
            break;
        }
    }

    return_if_valid_else!(ans, {
        check!(@nan, elliprd, [x, y, z]);
        case!(@any [x, y, z] == inf!(), T::zero());
        Err("elliprd: Failed to converge.")
    })
}

#[cfg(not(feature = "reduce-iteration"))]
const N_MAX_ITERATIONS: usize = 50;

#[cfg(feature = "reduce-iteration")]
const N_MAX_ITERATIONS: usize = 1;

#[cfg(not(feature = "reduce-iteration"))]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::compare_test_data_boost;

    fn _elliprd(inp: &[f64]) -> f64 {
        elliprd(inp[0], inp[1], inp[2]).unwrap()
    }

    #[test]
    fn test_elliprd() {
        compare_test_data_boost!("elliprd_data.txt", _elliprd, 4.8e-16);
    }

    #[test]
    fn test_elliprd_0xy() {
        compare_test_data_boost!("elliprd_0xy.txt", _elliprd, 5.9e-16);
    }

    #[test]
    fn test_elliprd_0yy() {
        compare_test_data_boost!("elliprd_0yy.txt", _elliprd, 2.6e-16);
    }
    #[test]
    fn test_elliprd_xxx() {
        compare_test_data_boost!("elliprd_xxx.txt", _elliprd, 2.3e-16);
    }

    #[test]
    fn test_elliprd_xxz() {
        compare_test_data_boost!("elliprd_xxz.txt", _elliprd, 7.9e-16);
    }

    #[test]
    fn test_elliprd_xyy() {
        compare_test_data_boost!("elliprd_xyy.txt", _elliprd, 3.7e-15);
    }

    #[test]
    fn test_elliprd_special_cases() {
        use std::f64::{INFINITY, NAN};
        // x < 0 or y < 0: should return Err
        assert!(elliprd(-1.0, 1.0, 1.0).is_err());
        assert!(elliprd(1.0, -1.0, 1.0).is_err());
        // z <= 0: should return Err
        assert!(elliprd(1.0, 1.0, 0.0).is_err());
        assert!(elliprd(1.0, 1.0, -1.0).is_err());
        // both x and y zero: should return Err
        assert!(elliprd(0.0, 0.0, 1.0).is_err());
        // NANs: should return Err
        assert!(elliprd(NAN, 1.0, 1.0).is_err());
        assert!(elliprd(1.0, NAN, 1.0).is_err());
        assert!(elliprd(1.0, 1.0, NAN).is_err());
        // Infs: should return zero
        assert_eq!(elliprd(INFINITY, 1.0, 1.0).unwrap(), 0.0);
        assert_eq!(elliprd(1.0, INFINITY, 1.0).unwrap(), 0.0);
        assert_eq!(elliprd(1.0, 1.0, INFINITY).unwrap(), 0.0);
    }
}

#[cfg(feature = "reduce-iteration")]
crate::test_force_unreachable! {
    assert!(elliprd(0.2, 0.5, 1e300).is_err());
}
