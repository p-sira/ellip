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
            "Assertion failed: expected = {:?}, got = {:?}, relative = {:?}, rtol = {:?}",
            expected, actual, relative, rtol
        )
    }
}
