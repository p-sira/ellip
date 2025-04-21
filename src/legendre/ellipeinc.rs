/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Boost Math and SciPy, see LICENSE in this directory.
 */

// Original header from Boost Math
//  Copyright (c) 2006 Xiaogang Zhang
//  Copyright (c) 2006 John Maddock
//  Copyright (c) 2024 Matt Borland
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  History:
//  XZ wrote the original of this file as part of the Google
//  Summer of Code 2006.  JM modified it to fit into the
//  Boost.Math conceptual framework better, and to ensure
//  that the code continues to work no matter how many digits
//  type T has.

use num_traits::Float;

use crate::{ellipe, elliprd, elliprf};

/// Computes [incomplete elliptic integral of the second kind](https://dlmf.nist.gov/19.2.E5).
/// ```text
///              φ
///             ⌠   _____________
/// E(φ, m)  =  │ \╱ 1 - m sin²θ  dθ
///             ⌡
///            0
/// ```
///
/// ## Parameters
/// - phi: amplitude angle (φ). φ ∈ ℝ.
/// - m: elliptic parameter. m ∈ ℝ.
///
/// The elliptic modulus (k) is also frequently used instead of the parameter (m), where k² = m.
///
/// ## Domain
/// - Returns error if m sin²φ > 1.
///
/// ## Graph
/// ![Incomplete Elliptic Integral of the Second Kind](https://github.com/p-sira/ellip/blob/main/figures/ellipeinc_plot.svg?raw=true)
///
/// [Interactive Plot](https://github.com/p-sira/ellip/blob/main/figures/ellipeinc_plot.html)
///
/// ![Incomplete Elliptic Integral of the Second Kind (3D Plot)](https://github.com/p-sira/ellip/blob/main/figures/ellipeinc_plot_3d.png?raw=true)
///
/// [Interactive 3D Plot](https://github.com/p-sira/ellip/blob/main/figures/ellipeinc_plot_3d.html)
///
/// # Related Functions
/// With c = csc²φ,
/// - [ellipeinc](crate::ellipeinc)(φ, m) = [elliprf](crate::elliprf)(c - 1, c - m, c) - m / 3 * [elliprd](crate::elliprd)(c - 1, c - m, c)
///
/// # Examples
/// ```
/// use ellip::{ellipeinc, util::assert_close};
/// use std::f64::consts::FRAC_PI_4;
///
/// assert_close(ellipeinc(FRAC_PI_4, 0.5).unwrap(), 0.7481865041776612, 1e-15);
/// ```
///
/// # References
/// - Maddock, John, Paul Bristow, Hubert Holin, and Xiaogang Zhang. “Boost Math Library: Special Functions - Elliptic Integrals.” Accessed April 17, 2025. <https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/ellint.html>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
/// - The MathWorks, Inc. “ellipticE.” Accessed April 21, 2025. <https://www.mathworks.com/help/symbolic/sym.elliptice.html>.
///
pub fn ellipeinc<T: Float>(phi: T, m: T) -> Result<T, &'static str> {
    if phi == zero!() {
        return Ok(zero!());
    }

    if m == zero!() {
        return Ok(phi);
    }

    let (phi, invert) = if phi < zero!() {
        (-phi, -one!())
    } else {
        (phi, one!())
    };

    if phi >= max_val!() {
        return Ok(invert * inf!());
    }

    if phi > one!() / epsilon!() {
        return Ok(invert * two!() * phi * ellipe(m)? / pi!());
    }

    if m == one!() {
        // For k = 1 ellipse actually turns to a line and every pi/2 in phi is exactly 1 in arc length
        // Periodicity though is in pi, curve follows sin(pi) for 0 <= phi <= pi/2 and then
        // 2 - sin(pi- phi) = 2 + sin(phi - pi) for pi/2 <= phi <= pi, so general form is:
        //
        // 2n + sin(phi - n * pi) ; |phi - n * pi| <= pi / 2
        let mm = (phi / pi!()).round();
        let remains = phi - mm * pi!();
        return Ok(invert * (two!() * mm + remains.sin()));
    }

    // Carlson's algorithm works only for |phi| <= pi/2,
    // use the integrand's periodicity to normalize phi
    //
    // Xiaogang's original code used a cast to long long here
    // but that fails if T has more digits than a long long,
    // so rewritten to use fmod instead:
    let mut rphi = phi % pi_2!();
    let mut mm = ((phi - rphi) / pi_2!()).round();
    let mut s = one!();
    if mm % two!() > half!() {
        mm = mm + one!();
        s = -one!();
        rphi = pi_2!() - rphi;
    }

    let mut result = if m > zero!() && rphi.powi(3) * m / six!() < epsilon!() * rphi.abs() {
        // See http://functions.wolfram.com/EllipticIntegrals/EllipticE2/06/01/03/0001/
        s * rphi
    } else {
        let s2p = rphi.sin() * rphi.sin();
        if m * s2p >= one!() {
            return Err("ellipeinc: m sin²φ must be smaller than one.");
        }
        let c2p = rphi.cos() * rphi.cos();
        let c = one!() / s2p;
        let cm1 = c2p / s2p;
        s * ((one!() - m) * elliprf(cm1, c - m, c)?
            + m * (one!() - m) * elliprd(cm1, c, c - m)? / three!()
            + m * (cm1 / (c * (c - m))).sqrt())
    };

    if mm != zero!() {
        result = result + mm * ellipe(m)?;
    }

    Ok(invert * result)
}

#[cfg(test)]
mod tests {
    use crate::compare_test_data;

    use super::*;

    fn ellipeinc_k<T: Float>(inp: &[T]) -> T {
        ellipeinc(inp[0], inp[1] * inp[1]).unwrap()
    }

    #[test]
    fn test_ellipeinc() {
        compare_test_data!(
            "./tests/data/boost/ellipeinc_data.txt",
            ellipeinc_k::<f64>,
            5e-16
        );
    }
}
