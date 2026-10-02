/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is translated from SciPy C++ implementation to Rust.
 */

/* Translated into C++ by SciPy developers in 2024. */

/*
 * Cephes Math Library Release 2.1:  December, 1988
 * Copyright 1984, 1987, 1988 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

/* Sources:
 * [1] Holin et. al., "Polynomial and Rational Function Evaluation",
 *     https://www.boost.org/doc/libs/1_61_0/libs/math/doc/html/math_toolkit/roots/rational.html
 */

use num_traits::Float;

/// Evaluate polynomial with coefficients in reverse order (C_0 + C_1 x + C_2 x^2 + ...)
#[inline]
pub(crate) fn polyeval<T: Float>(x: T, coeff: &[T]) -> T {
    if coeff.len() == 12 {
        let x2 = x * x;
        let x4 = x2 * x2;
        let x8 = x4 * x4;

        let p0 = coeff[0] + coeff[1] * x;
        let p1 = coeff[2] + coeff[3] * x;
        let p2 = coeff[4] + coeff[5] * x;
        let p3 = coeff[6] + coeff[7] * x;
        let p4 = coeff[8] + coeff[9] * x;
        let p5 = coeff[10] + coeff[11] * x;

        let q0 = p0 + p1 * x2;
        let q1 = p2 + p3 * x2;
        let q2 = p4 + p5 * x2;

        let r0 = q0 + q1 * x4;
        return r0 + q2 * x8;
    }

    if coeff.len() == 6 {
        let x2 = x * x;
        let x4 = x2 * x2;

        let p0 = coeff[0] + coeff[1] * x;
        let p1 = coeff[2] + coeff[3] * x;
        let p2 = coeff[4] + coeff[5] * x;

        let q0 = p0 + p1 * x2;
        return q0 + p2 * x4;
    }

    let mut ans = T::zero();
    coeff.iter().rev().for_each(|&k| ans = ans * x + k);
    ans
}
