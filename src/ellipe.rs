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

use crate::polyeval;

/// Compute complete elliptic integral of the second kind.
///
/// ```text
///           π/2
///          ⌠     _______________
/// E(m)  =  │   \╱ 1 - m sin²(t)  dt
///          ⌡
///         0
/// where m ≤ 1
/// ```
///
pub fn ellipe(m: f64) -> Result<f64, &'static str> {
    if m > 1.0 {
        return Err("ellipe: m must satisfy: m ≤ 1.");
    }

    if m == 1.0 {
        return Ok(1.0);
    }

    let mut m = m;
    let mut k = 1.0;
    while m < 0.0 {
        k *= (1.0 - m).sqrt();
        m /= m - 1.0;
    }

    let x: f64 = 1.0 - m;

    let p_val = polyeval(x, &ELLPE_P, 10);
    let log_x = x.ln();

    Ok(k * (p_val - log_x * (x * polyeval(x, &ELLPE_Q, 9))))
}

const ELLPE_P: [f64; 11] = [
    1.535_525_773_010_133E-4,
    2.508_884_921_636_020_4E-3,
    8.687_868_165_658_896E-3,
    1.073_509_490_560_761_9E-2,
    7.773_954_925_167_871E-3,
    7.583_952_894_135_147E-3,
    1.156_884_368_105_741_2E-2,
    2.183_179_960_155_572_4E-2,
    5.680_519_456_178_606E-2,
    4.431_471_805_609_908_4E-1,
    1.0,
];

const ELLPE_Q: [f64; 10] = [
    3.279_548_985_764_858_5E-5,
    1.009_627_926_793_567_2E-3,
    6.506_094_899_769_275E-3,
    1.688_621_639_933_113_3E-2,
    2.617_697_424_544_936_4E-2,
    3.348_339_048_882_249E-2,
    4.271_809_265_189_315E-2,
    5.859_366_344_711_01E-2,
    9.374_999_971_976_443E-2,
    2.499_999_999_998_883E-1,
];

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
