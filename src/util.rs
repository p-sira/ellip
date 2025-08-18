/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

//! Utility functions like assert_close.

use num_traits::Float;
use std::fmt::Debug;

/// Assert that the actual value is within the relative tolerance of the expected value.
///
/// Panics if the assertion is failed.
pub fn assert_close<T: Float + Debug>(actual: T, expected: T, rtol: T) {
    let relative = (actual - expected).abs() / expected;
    if relative > rtol {
        panic!(
            "Assertion failed: expected = {expected:?}, got = {actual:?}, relative = {relative:?}, rtol = {rtol:?}",
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
}
