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

use num_traits::Float;

use crate::{
    crate_util::{case, check, let_mut},
    StrErr,
};

/// Computes RC ([degenerate symmetric elliptic integral of RF](https://dlmf.nist.gov/19.16.E6)).
/// ```text
///                  ∞                  
///              1  ⌠        dt        
/// RC(x, y)  =  ─  ⎮ ─────────────────
///              2  ⎮             _____
///                 ⌡ (t + y) ⋅ ╲╱t + x
///                0                  
/// ```
///
/// ## Parameters
/// - x ∈ ℝ, x ≥ 0
/// - y ∈ ℝ, y ≠ 0
///
/// ## Domain
/// - Returns error if x < 0 or y = 0.
/// - Returns the Cauchy principal value if y < 0.
///
/// ## Graph
/// ![Degenerate Symmetric Elliptic Integral of RF (RC)](https://github.com/p-sira/ellip/blob/main/figures/elliprc_plot.svg?raw=true)
///
/// [Interactive Plot](https://github.com/p-sira/ellip/blob/main/figures/elliprc_plot.html)
///
/// ## Special Cases
/// - RC(x, x) = 1/sqrt(x)
/// - RC(0, y) = π/(2*sqrt(y))
/// - RC(x, y) = atan(sqrt(y-x)/x) / sqrt(y-x) for y > x
/// - RC(x, y) = ln(sqrt(x) + sqrt(x-y)) / sqrt(x-y) for y < x
/// - RC(x, y) = 0 for x = ∞ or y = ∞
///
/// ## Notes
/// RC is a degenerate case of the RF. It is an elementary function rather than an elliptic integral.
///
/// # Related Functions
/// - [elliprc](crate::elliprc)(x, y) = [elliprf](crate::elliprf)(x, y, y)
///
/// # Examples
/// ```
/// use ellip::{elliprc, util::assert_close};
///
/// assert_close(elliprc(1.0, 0.5).unwrap(), 1.2464504802804608, 1e-15);
/// ```
///
/// # References
/// - Maddock, John, Paul Bristow, Hubert Holin, and Xiaogang Zhang. “Boost Math Library: Special Functions - Elliptic Integrals.” Accessed April 17, 2025. <https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/ellint.html>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
/// - The SciPy Community. “SciPy: Special Functions - Elliprc.” Accessed April 17, 2025. <https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.elliprc.html>.
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn elliprc<T: Float>(x: T, y: T) -> Result<T, StrErr> {
    if x < 0.0 {
        return Err("elliprc: x must be non-negative.");
    }

    if y == 0.0 {
        return Err("elliprc: y must be non-zero.");
    }

    let ans = elliprc_unchecked(x, y);
    if ans.is_finite() {
        return Ok(ans);
    }
    check!(@nan, elliprc, [x, y]);
    case!(@any [x, y] == inf!(), T::zero());
    Err("elliprc: Unexpected error.")
}

/// Unsafe version of [elliprc](crate::elliprc).
///
/// Undefined behavior with invalid arguments and edge cases.
/// # Known Invalid Cases
/// - x < 0
/// - y = 0
/// - x = ∞ or y = ∞
#[inline]
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn elliprc_unchecked<T: Float>(x: T, y: T) -> T {
    let_mut!(x, y);
    let mut prefix = 1.0;
    // for y < 0, the integral is singular, return Cauchy principal value
    if y < 0.0 {
        prefix = (x / (x - y)).sqrt();
        x = x - y;
        y = -y;
    }

    if x == 0.0 {
        return prefix * pi_2!() / y.sqrt();
    }

    if x == y {
        return prefix / x.sqrt();
    }

    if y > x {
        return prefix * ((y - x) / x).sqrt().atan() / (y - x).sqrt();
    }

    #[cfg(feature = "reduce-iteration")]
    return nan!();

    if y / x > 0.5 {
        let arg = ((x - y) / x).sqrt();
        prefix * ((arg).ln_1p() - (-arg).ln_1p()) / (2.0 * (x - y).sqrt())
    } else {
        prefix * ((x.sqrt() + (x - y).sqrt()) / y.sqrt()).ln() / (x - y).sqrt()
    }
}

#[cfg(not(feature = "reduce-iteration"))]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::compare_test_data_boost;

    fn _elliprc(inp: &[f64]) -> f64 {
        elliprc(inp[0], inp[1]).unwrap()
    }

    #[test]
    fn test_elliprc() {
        compare_test_data_boost!("elliprc_data.txt", _elliprc, f64::EPSILON);
    }

    #[test]
    fn test_elliprc_special_cases() {
        use std::f64::{consts::PI, INFINITY, NAN};
        // x < 0: should return Err
        assert_eq!(elliprc(-1.0, 1.0), Err("elliprc: x must be non-negative."));
        // y == 0: should return Err
        assert_eq!(elliprc(1.0, 0.0), Err("elliprc: y must be non-zero."));
        // RC(x, x) = 1/sqrt(x)
        assert_eq!(elliprc(1.0, 1.0).unwrap(), 1.0);
        assert_eq!(elliprc(4.0, 4.0).unwrap(), 0.5);
        // RC(0, y) = π/(2*sqrt(y))
        assert_eq!(elliprc(0.0, 1.0).unwrap(), PI / 2.0);
        assert_eq!(elliprc(0.0, 4.0).unwrap(), PI / 4.0);
        // RC(x, y) = atan(sqrt(y-x)/x) / sqrt(y-x) for y > x
        assert_eq!(elliprc(1.0, 4.0).unwrap(), (3.0.sqrt().atan() / 3.0.sqrt()));
        assert_eq!(elliprc(1.0, 9.0).unwrap(), (8.0.sqrt().atan() / 8.0.sqrt()));
        // RC(x, y) = ln(sqrt(x) + sqrt(x-y)) / sqrt(x-y) for y < x
        assert_eq!(
            elliprc(4.0, 1.0).unwrap(),
            ((2.0 + 3.0.sqrt()).ln() / 3.0.sqrt())
        );
        assert_eq!(
            elliprc(9.0, 1.0).unwrap(),
            ((3.0 + 8.0.sqrt()).ln() / 8.0.sqrt())
        );
        // NANs: should return Err
        assert_eq!(elliprc(NAN, 1.0), Err("elliprc: Arguments cannot be NAN."));
        assert_eq!(elliprc(1.0, NAN), Err("elliprc: Arguments cannot be NAN."));
        // Infs: should return 0
        assert_eq!(elliprc(INFINITY, 1.0).unwrap(), 0.0);
        assert_eq!(elliprc(1.0, INFINITY).unwrap(), 0.0);
    }
}

#[cfg(feature = "reduce-iteration")]
crate::test_force_unreachable! {
    assert_eq!(elliprc(2.0, 1.0), Err("elliprc: Unexpected error."));
}
