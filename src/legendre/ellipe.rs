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
    if m > one!() {
        return Err("ellipe: m must be less than 1.");
    }

    // Special cases: https://dlmf.nist.gov/19.6.E1
    if m == one!() {
        return Ok(one!());
    }

    if m == zero!() {
        return Ok(pi_2!());
    }

    // Note: this function allows both negative m and positive m less than 1.
    // Negative m: Abramowitz & Stegun, 1972
    let mut m = m;
    let mut c = one!();
    while m < zero!() {
        c = c * (one!() - m).sqrt();
        m = m / (m - one!());
    }

    let x = one!() - m;

    Ok(c * polyeval(x, &ellpe_p()) - x.ln() * (x * polyeval(x, &ellpe_q())))
}

#[inline]
fn ellpe_p<T: Float>() -> [T; 11] {
    [
        num!(1.53552577301013293365E-4),
        num!(2.50888492163602060990E-3),
        num!(8.68786816565889628429E-3),
        num!(1.07350949056076193403E-2),
        num!(7.77395492516787092951E-3),
        num!(7.58395289413514708519E-3),
        num!(1.15688436810574127319E-2),
        num!(2.18317996015557253103E-2),
        num!(5.68051945617860553470E-2),
        num!(4.43147180560990850618E-1),
        num!(1.00000000000000000299E0),
    ]
}

#[inline]
fn ellpe_q<T: Float>() -> [T; 10] {
    [
        num!(3.27954898576485872656E-5),
        num!(1.00962792679356715133E-3),
        num!(6.50609489976927491433E-3),
        num!(1.68862163993311317300E-2),
        num!(2.61769742454493659583E-2),
        num!(3.34833904888224918614E-2),
        num!(4.27180926518931511717E-2),
        num!(5.85936634471101055642E-2),
        num!(9.37499997197644278445E-2),
        num!(2.49999999999888314361E-1),
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
        compare_test_data!("./tests/data/boost/ellint_e_data.txt", ellipe_k, RTOL);
    }
}
