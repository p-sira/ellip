/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is translated from SciPy C++ implementation to Rust.
 * Modified to use iteration instead of recursion.
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

use crate::polyeval;

/// Compute [complete elliptic integral of the first kind](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipk.html).
///
/// ```text
///           π/2
///          ⌠          dt
/// F(m)  =  │  ___________________
///          │     _______________
///          ⌡   \╱ 1 - m sin²(t)
///         0
/// where m ≤ 1
/// ```
///
/// Note that some mathematical references use the parameter k for the function,
/// where k² = m.
pub fn ellipk(m: f64) -> Result<f64, &'static str> {
    if m > 1.0 {
        return Err("ellipk: m must be less than 1.");
    }

    let mut x = 1.0 - m;

    if x.is_infinite() {
        return Ok(0.0);
    }

    let mut k = 1.0;
    while x > 1.0 {
        k /= x.sqrt();
        x = 1.0 / x;
    }

    if x > f64::EPSILON {
        return Ok(k * (polyeval(x, &ELLPK_P, 10) - x.ln() * polyeval(x, &ELLPK_Q, 10)));
    }

    if x == 0.0 {
        return Ok(f64::INFINITY);
    }

    Ok(k * (4.0_f64.ln() - 0.5 * x.ln()))
}

const ELLPK_P: [f64; 11] = [
    1.379_828_646_062_732_5E-4,
    2.280_257_240_058_756E-3,
    7.974_040_132_204_152E-3,
    9.858_213_790_212_26E-3,
    6.874_896_874_499_499E-3,
    6.189_010_336_376_876E-3,
    8.790_782_739_527_438E-3,
    1.493_804_489_168_052_6E-2,
    3.088_514_652_467_12E-2,
    9.657_359_028_116_902E-2,
    1.386_294_361_119_890_6,
];

const ELLPK_Q: [f64; 11] = [
    2.940_789_550_485_985E-5,
    9.141_847_238_659_173E-4,
    5.940_583_037_531_678E-3,
    1.548_505_166_497_624E-2,
    2.390_896_027_159_248_8E-2,
    3.012_047_152_276_040_4E-2,
    3.737_743_141_738_232_6E-2,
    4.882_803_475_709_983E-2,
    7.031_249_969_639_575E-2,
    1.249_999_999_998_708_3E-1,
    5E-1,
];

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{compare_test_data, test_util::RTOL};

    fn ellipk_k(k: &[f64]) -> f64 {
        ellipk(k[0] * k[0]).unwrap()
    }

    #[test]
    fn test_ellipk() {
        compare_test_data!("./tests/data/boost/ellint_k_data.txt", ellipk_k, RTOL);
    }
}
