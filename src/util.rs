/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

//! Utility functions like assert_close.

use num_traits::Float;

/// Assert that the actual value is within the relative tolerance of the expected value.
///
/// Panics if the assertion is failed.
pub fn assert_close<T: Float>(actual: T, expected: T, rtol: T) {
    let relative = (actual - expected).abs() / expected.abs();
    let valid_tolerance = rtol.is_finite() && rtol >= T::zero();
    let close =
        actual == expected || (actual.is_finite() && expected.is_finite() && relative <= rtol);
    if !valid_tolerance || !close {
        panic!(
            "Assertion failed: expected = {}, got = {}, relative = {}, rtol = {}",
            expected.to_f64().unwrap(),
            actual.to_f64().unwrap(),
            relative.to_f64().unwrap(),
            rtol.to_f64().unwrap()
        )
    }
}

#[cfg(not(feature = "test_force_fail"))]
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic]
    fn test_assert_close_panic() {
        assert_close(1.0, 2.0, 1e-6);
    }

    #[test]
    fn test_assert_close_success() {
        assert_close(1.0, 1.0 + 1e-6, 1e-6);
    }

    // Regression for audit finding A15: https://github.com/p-sira/ellip/pull/129
    #[test]
    fn test_reject_invalid_comparisons() {
        for (actual, expected, tol) in [
            (100.0, -1.0, 1e-15),
            (f64::NAN, 1.0, 1e-15),
            (1.0, f64::NAN, 1e-15),
            (1.0, 0.0, 1e-15),
            (f64::INFINITY, 1.0, 1e-15),
            (1.0, f64::INFINITY, 1e-15),
            (f64::INFINITY, f64::NEG_INFINITY, 1e-15),
            (1.0, 1.0, f64::NAN),
            (1.0, 1.0, -1.0),
        ] {
            assert!(
                std::panic::catch_unwind(|| crate::util::assert_close(actual, expected, tol))
                    .is_err(),
                "accepted {actual}, {expected}, {tol}"
            );
        }
    }
    #[test]
    fn test_accept_exact_and_signed_values() {
        for x in [-1.0, 0.0, -0.0, 1.0, f64::INFINITY, f64::NEG_INFINITY] {
            crate::util::assert_close(x, x, 0.0);
        }
        crate::util::assert_close(-1.0 - 1e-8, -1.0, 2e-8);
    }
}
