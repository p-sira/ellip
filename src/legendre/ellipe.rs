/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is translated from SciPy C++ implementation to Rust.
 * Modified to use iteration instead of recursion.
 */

/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     ellpe.c
 *
 *     Complete elliptic integral of the second kind
 *
 *
 *
 * SYNOPSIS:
 *
 * double m, y, ellpe();
 *
 * y = ellpe( m );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *            pi/2
 *             -
 *            | |                 2
 * E(m)  =    |    sqrt( 1 - m sin t ) dt
 *          | |
 *           -
 *            0
 *
 * Where m = 1 - m1, using the approximation
 *
 *      P(x)  -  x log x Q(x).
 *
 * Though there are no singularities, the argument m1 is used
 * internally rather than m for compatibility with ellpk().
 *
 * E(1) = 1; E(0) = pi/2.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0, 1       10000       2.1e-16     7.3e-17
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * ellpe domain      x<0, x>1            0.0
 *
 */

/*                                                     ellpe.c         */

/* Elliptic integral of second kind */

/*
 * Cephes Math Library, Release 2.1:  February, 1989
 * Copyright 1984, 1987, 1989 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 *
 * Feb, 2002:  altered by Travis Oliphant
 * so that it is called with argument m
 * (which gets immediately converted to m1 = 1-m)
 */

use std::f64::consts::FRAC_PI_2;

use num_traits::Float;

use crate::polyeval;

/// Compute [complete elliptic integral of the second kind](https://dlmf.nist.gov/19.2.E8).
///
/// ```text
///           π/2
///          ⌠     ___________
/// E(m)  =  │ \╱ 1 - m sin²θ  dθ
///          ⌡
///         0
/// where m ≤ 1
/// ```
///
/// Note that some mathematical references use the parameter k for the function,
/// where k² = m.
pub fn ellipe<T: Float>(m: T) -> Result<T, &'static str> {
    if m > T::one() {
        return Err("ellipe: m must be less than 1.");
    }

    // Special cases: https://dlmf.nist.gov/19.6.E1
    if m == T::one() {
        return Ok(T::one());
    }

    if m == T::zero() {
        return Ok(T::from(FRAC_PI_2).unwrap());
    }

    // Note: this function allows both negative m and positive m less than 1.
    Ok(_ellipe_neg(m))
}

/// Unchecked version of [ellipe].
///
/// Domain: 0 < m ≤ 1
pub fn _ellipe<T: Float>(m: T) -> T {
    let x = T::one() - m;

    polyeval(x, &ellpe_p()) - x.ln() * (x * polyeval(x, &ellpe_q()))
}

/// [_ellipe] but allows m < 0. An unchecked version of [ellipe].
///
/// Domain: m ≤ 1
pub fn _ellipe_neg<T: Float>(m: T) -> T {
    // Negative m: Abramowitz & Stegun, 1972
    let mut m = m;
    let mut c = T::one();
    while m < T::zero() {
        c = c * (T::one() - m).sqrt();
        m = m / (m - T::one());
    }

    let x = T::one() - m;

    c * polyeval(x, &ellpe_p()) - x.ln() * (x * polyeval(x, &ellpe_q()))
}

#[inline]
fn ellpe_p<T: Float>() -> [T; 11] {
    [
        T::from(1.53552577301013293365E-4).unwrap(),
        T::from(2.50888492163602060990E-3).unwrap(),
        T::from(8.68786816565889628429E-3).unwrap(),
        T::from(1.07350949056076193403E-2).unwrap(),
        T::from(7.77395492516787092951E-3).unwrap(),
        T::from(7.58395289413514708519E-3).unwrap(),
        T::from(1.15688436810574127319E-2).unwrap(),
        T::from(2.18317996015557253103E-2).unwrap(),
        T::from(5.68051945617860553470E-2).unwrap(),
        T::from(4.43147180560990850618E-1).unwrap(),
        T::from(1.00000000000000000299E0).unwrap(),
    ]
}

#[inline]
fn ellpe_q<T: Float>() -> [T; 10] {
    [
        T::from(3.27954898576485872656E-5).unwrap(),
        T::from(1.00962792679356715133E-3).unwrap(),
        T::from(6.50609489976927491433E-3).unwrap(),
        T::from(1.68862163993311317300E-2).unwrap(),
        T::from(2.61769742454493659583E-2).unwrap(),
        T::from(3.34833904888224918614E-2).unwrap(),
        T::from(4.27180926518931511717E-2).unwrap(),
        T::from(5.85936634471101055642E-2).unwrap(),
        T::from(9.37499997197644278445E-2).unwrap(),
        T::from(2.49999999999888314361E-1).unwrap(),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{compare_test_data, test_util::RTOL};

    fn ellipe_k(k: &[f64]) -> f64 {
        ellipe(k[0] * k[0]).unwrap()
    }

    #[test]
    fn test_ellipe() {
        compare_test_data!("./tests/data/boost/ellipe_data.txt", ellipe_k, RTOL);
    }
}
