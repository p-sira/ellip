/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is translated from SciPy C++ implementation to Rust.
 * Modified to use iteration instead of recursion and some special case handling.
 */

/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     ellpk.c
 *
 *     Complete elliptic integral of the first kind
 *
 *
 *
 * SYNOPSIS:
 *
 * double m1, y, ellpk();
 *
 * y = ellpk( m1 );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *
 *            pi/2
 *             -
 *            | |
 *            |           dt
 * K(m)  =    |    ------------------
 *            |                   2
 *          | |    sqrt( 1 - m sin t )
 *           -
 *            0
 *
 * where m = 1 - m1, using the approximation
 *
 *     P(x)  -  log x Q(x).
 *
 * The argument m1 is used internally rather than m so that the logarithmic
 * singularity at m = 1 will be shifted to the origin; this
 * preserves maximum accuracy.
 *
 * K(0) = pi/2.
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,1        30000       2.5e-16     6.8e-17
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * ellpk domain       x<0, x>1           0.0
 *
 */

/*                                                     ellpk.c */

/*
 * Cephes Math Library, Release 2.0:  April, 1987
 * Copyright 1984, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

use num_traits::Float;

use crate::polyeval;

/// Compute [complete elliptic integral of the first kind](https://dlmf.nist.gov/19.2.E8).
/// ```text
///           π/2
///          ⌠          dθ
/// K(m)  =  │  _________________
///          │     _____________
///          ⌡   \╱ 1 - m sin²θ
///         0
/// where m ≤ 1
/// ```
///
/// Note that some mathematical references use the parameter k for the function,
/// where k² = m.
///
/// # Examples
/// ```
/// use ellip::{ellipk, util::assert_close};
///
/// assert_close(ellipk(0.5).unwrap(), 1.8540746773013719, 1e-15);
/// ```
pub fn ellipk<T: Float>(m: T) -> Result<T, &'static str> {
    if m > one!() {
        return Err("ellipk: m must be less than 1.");
    }

    // Special cases: https://dlmf.nist.gov/19.6.E1
    if m == one!() {
        return Ok(inf!());
    }

    if m == zero!() {
        return Ok(pi_2!());
    }

    let x: T = one!() - m;

    // Asymptotic approximation when m -> 1
    if x < epsilon!() {
        return Ok(four!().ln() - half!() * x.ln());
    }

    Ok(polyeval(x, &ellpk_p::<T>()) - x.ln() * polyeval(x, &ellpk_q::<T>()))
}

#[inline]
fn ellpk_p<T: Float>() -> [T; 11] {
    [
        num!(1.37982864606273237150E-4),
        num!(2.28025724005875567385E-3),
        num!(7.97404013220415179367E-3),
        num!(9.85821379021226008714E-3),
        num!(6.87489687449949877925E-3),
        num!(6.18901033637687613229E-3),
        num!(8.79078273952743772254E-3),
        num!(1.49380448916805252718E-2),
        num!(3.08851465246711995998E-2),
        num!(9.65735902811690126535E-2),
        num!(1.38629436111989062502E0),
    ]
}

#[inline]
fn ellpk_q<T: Float>() -> [T; 11] {
    [
        num!(2.94078955048598507511E-5),
        num!(9.14184723865917226571E-4),
        num!(5.94058303753167793257E-3),
        num!(1.54850516649762399335E-2),
        num!(2.39089602715924892727E-2),
        num!(3.01204715227604046988E-2),
        num!(3.73774314173823228969E-2),
        num!(4.88280347570998239232E-2),
        num!(7.03124996963957469739E-2),
        num!(1.24999999999870820058E-1),
        num!(4.99999999999999999821E-1),
    ]
}

#[cfg(test)]
mod tests {
    use core::f64;

    use super::*;
    use crate::compare_test_data;

    fn ellipk_k(k: &[f64]) -> f64 {
        ellipk(k[0] * k[0]).unwrap()
    }

    #[test]
    fn test_ellipk() {
        compare_test_data!("./tests/data/boost/ellipk_data.txt", ellipk_k, f64::EPSILON);
    }
}
