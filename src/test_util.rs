/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

#[macro_export]
macro_rules! assert_close {
    ($expected: expr, $actual: expr, $rtol: expr) => {
        $crate::util::assert_close($actual, $expected, $rtol);
    };
}

#[cfg(feature = "test_force_fail")]
#[macro_export]
macro_rules! test_force_unreachable {
    ($($inner:tt)*) => {
        #[cfg(test)]
        mod tests {
            use super::*;

            #[test]
            fn force_unreachable() {
                $($inner)*
            }
        }
    };
}

#[cfg(all(test, not(feature = "test_force_fail")))]
mod assertion_regression_tests {
    // Regression for audit finding A15: https://github.com/p-sira/ellip/pull/129
    #[test]
    fn macro_rejects_negative_reference_errors_and_nan() {
        assert!(std::panic::catch_unwind(|| {
            crate::assert_close!(-1.0_f64, 100.0_f64, 1e-15);
        })
        .is_err());
        assert!(std::panic::catch_unwind(|| {
            crate::assert_close!(f64::NAN, 1.0_f64, 1e-15);
        })
        .is_err());
        crate::assert_close!(-1.0_f64, -1.0_f64, 0.0);
        crate::assert_close!(0.0_f64, 0.0_f64, 0.0);
    }
}
