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
    let mut ans = T::zero();
    coeff.iter().rev().for_each(|&k| ans = ans * x + k);
    ans
}
