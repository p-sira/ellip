/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use num_traits::Float;
use std::fmt::Debug;

pub fn assert_close<T: Float + Debug>(expected: T, actual: T, rtol: T) {
    let relative = (actual - expected).abs() / expected;
    if relative > rtol {
        panic!(
            "Assertion failed: expected = {:?}, got = {:?}, relative = {:?}, rtol = {:?}",
            expected, actual, relative, rtol
        )
    }
}
