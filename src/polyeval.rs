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

    if coeff.len() == 9 {
        let x2 = x * x;
        let x4 = x2 * x2;
        let x8 = x4 * x4;

        let p0 = coeff[0] + coeff[1] * x;
        let p1 = coeff[2] + coeff[3] * x;
        let p2 = coeff[4] + coeff[5] * x;
        let p3 = coeff[6] + coeff[7] * x;

        let q0 = p0 + p1 * x2;
        let q1 = p2 + p3 * x2;

        let r0 = q0 + q1 * x4;
        return r0 + coeff[8] * x8;
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

#[cfg(test)]
mod tests {
    use super::*;

    fn horner_eval<T: Float>(x: T, coeff: &[T]) -> T {
        let mut ans = T::zero();
        coeff.iter().rev().for_each(|&k| ans = ans * x + k);
        ans
    }

    #[test]
    fn test_polyeval_degree_12() {
        let coeff = [
            1.5, -2.0, 0.5, 3.2, -1.1, 0.8, -0.4, 1.2, -0.9, 0.3, -0.1, 0.05,
        ];
        for &x in &[0.0, 0.5, -0.5, 1.0, -1.0, 0.123, -0.876] {
            let actual = polyeval(x, &coeff);
            let expected = horner_eval(x, &coeff);
            assert!((actual - expected).abs() <= 1e-14 * expected.abs().max(1.0));
        }
    }

    #[test]
    fn test_polyeval_degree_8() {
        let coeff = [1.0, -0.33, 0.2, -0.14, 0.11, -0.09, 0.07, -0.06, 0.05];
        for &x in &[0.0, 0.5, -0.5, 1.0, -1.0, 0.01, -0.01, 0.35] {
            let actual = polyeval(x, &coeff);
            let expected = horner_eval(x, &coeff);
            assert!((actual - expected).abs() <= 1e-15 * expected.abs().max(1.0));
        }
    }

    #[test]
    fn test_polyeval_degree_6() {
        let coeff = [1.2, -0.8, 2.5, -1.4, 0.6, -0.2];
        for &x in &[0.0, 0.5, -0.5, 1.0, -1.0, 0.42, -0.73] {
            let actual = polyeval(x, &coeff);
            let expected = horner_eval(x, &coeff);
            assert!((actual - expected).abs() <= 1e-15 * expected.abs().max(1.0));
        }
    }

    #[test]
    fn test_polyeval_fallback_lengths() {
        // Empty
        assert_eq!(polyeval(2.0, &[]), 0.0);

        // Constant (len 1)
        assert_eq!(polyeval(2.0, &[3.5]), 3.5);

        // Linear (len 2)
        assert_eq!(polyeval(2.0, &[1.0, 2.0]), 5.0);

        // Quadratic (len 3)
        assert_eq!(polyeval(2.0, &[1.0, 2.0, 3.0]), 17.0);

        // Degree 4 (len 5)
        let coeff5 = [1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(polyeval(2.0, &coeff5), horner_eval(2.0, &coeff5));

        // Degree 6 (len 7)
        let coeff7 = [1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 4.0];
        assert_eq!(polyeval(1.5, &coeff7), horner_eval(1.5, &coeff7));

        // Degree 12 (len 13)
        let coeff13 = [1.0; 13];
        assert_eq!(polyeval(0.5, &coeff13), horner_eval(0.5, &coeff13));
    }

    #[test]
    fn test_polyeval_f32() {
        let coeff = [1.0f32, 2.0, 3.0, 4.0, 5.0, 6.0];
        assert_eq!(polyeval(0.5f32, &coeff), horner_eval(0.5f32, &coeff));

        let coeff12 = [0.1f32; 12];
        let diff = (polyeval(0.5f32, &coeff12) - horner_eval(0.5f32, &coeff12)).abs();
        assert!(diff <= 1e-6);
    }
}

