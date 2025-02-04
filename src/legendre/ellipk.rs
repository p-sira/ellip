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

use std::f64::consts::FRAC_PI_2;

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
pub fn ellipk<T: Float>(m: T) -> Result<T, &'static str> {
    if m > T::one() {
        return Err("ellipk: m must be less than 1.");
    }

    // Special cases: https://dlmf.nist.gov/19.6.E1
    if m == T::one() {
        return Ok(T::infinity());
    }

    if m == T::zero() {
        return Ok(T::from(FRAC_PI_2).unwrap());
    }

    let x: T = T::one() - m;

    // Asymptotic approximation when m -> 1
    if x < T::epsilon() {
        return Ok(T::from(4.0).unwrap().ln() - T::from(0.5).unwrap() * x.ln());
    }

    Ok(_ellipk(m))
}

/// Unchecked version of [ellipk].
///
/// Domain: m ≤ 1
pub fn _ellipk<T: Float>(m: T) -> T {
    let x = T::one() - m;
    polyeval(x, &ellpk_p::<T>()) - x.ln() * polyeval(x, &ellpk_q::<T>())
}

fn ellpk_p<T: Float>() -> [T; 11] {
    [
        T::from(1.37982864606273237150E-4).unwrap(),
        T::from(2.28025724005875567385E-3).unwrap(),
        T::from(7.97404013220415179367E-3).unwrap(),
        T::from(9.85821379021226008714E-3).unwrap(),
        T::from(6.87489687449949877925E-3).unwrap(),
        T::from(6.18901033637687613229E-3).unwrap(),
        T::from(8.79078273952743772254E-3).unwrap(),
        T::from(1.49380448916805252718E-2).unwrap(),
        T::from(3.08851465246711995998E-2).unwrap(),
        T::from(9.65735902811690126535E-2).unwrap(),
        T::from(1.38629436111989062502E0).unwrap(),
    ]
}

fn ellpk_q<T: Float>() -> [T; 11] {
    [
        T::from(2.94078955048598507511E-5).unwrap(),
        T::from(9.14184723865917226571E-4).unwrap(),
        T::from(5.94058303753167793257E-3).unwrap(),
        T::from(1.54850516649762399335E-2).unwrap(),
        T::from(2.39089602715924892727E-2).unwrap(),
        T::from(3.01204715227604046988E-2).unwrap(),
        T::from(3.73774314173823228969E-2).unwrap(),
        T::from(4.88280347570998239232E-2).unwrap(),
        T::from(7.03124996963957469739E-2).unwrap(),
        T::from(1.24999999999870820058E-1).unwrap(),
        T::from(4.99999999999999999821E-1).unwrap(),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{compare_test_data, test_util::RTOL};

    fn ellipk_k(k: &[f64]) -> f64 {
        ellipk(k[0] * k[0]).unwrap()
    }

    #[test]
    fn test_ellipk() {
        compare_test_data!("./tests/data/boost/ellipk_data.txt", ellipk_k, RTOL);
    }
}
