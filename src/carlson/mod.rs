/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

//! Elliptic integral functions in Carlson's form.

mod elliprc;
mod elliprd;
mod elliprf;
mod elliprg;
mod elliprj;

pub use elliprc::elliprc;
pub use elliprd::elliprd;
pub use elliprf::elliprf;
pub use elliprg::elliprg;
pub use elliprj::elliprj;

#[allow(unused_imports)]
#[cfg(not(feature = "unstable"))]
pub(crate) use {
    elliprc::elliprc_unchecked, elliprd::elliprd_unchecked, elliprf::elliprf_unchecked,
    elliprg::elliprf_rd_unchecked, elliprg::elliprg_unchecked, elliprj::elliprj_unchecked,
};

#[cfg(feature = "unstable")]
pub use {
    elliprc::elliprc_unchecked, elliprd::elliprd_unchecked, elliprf::elliprf_unchecked,
    elliprg::elliprf_rd_unchecked, elliprg::elliprg_unchecked, elliprj::elliprj_unchecked,
};

#[cfg(all(test, not(feature = "test_force_fail")))]
mod powf_regression_tests {
    // Paired numerical regression for the power replacements used by RD and
    // the combined RF/RD evaluation in RG. Include values on either side of
    // powers of two and extremes where (a * sqrt(a)) would overflow.
    #[test]
    fn inverse_three_halves_f64_matches_powf() {
        for exponent in (-680..=680).step_by(17) {
            let center = 2.0_f64.powi(exponent);
            for a in [
                f64::from_bits(center.to_bits() - 1),
                center,
                f64::from_bits(center.to_bits() + 1),
                1.37 * center,
            ] {
                let baseline = a.powf(-1.5);
                let optimized = a.sqrt().recip() / a;
                assert!(baseline.is_finite() && optimized.is_finite());
                assert!(
                    (optimized - baseline).abs() <= 8.0 * f64::EPSILON * baseline.abs(),
                    "a={a}, powf={baseline}, sqrt={optimized}"
                );
            }
        }
    }

    #[test]
    fn inverse_three_halves_f32_matches_powf() {
        for exponent in (-84..=84).step_by(3) {
            let center = 2.0_f32.powi(exponent);
            for a in [
                f32::from_bits(center.to_bits() - 1),
                center,
                f32::from_bits(center.to_bits() + 1),
                1.37 * center,
            ] {
                let baseline = a.powf(-1.5);
                let optimized = a.sqrt().recip() / a;
                assert!(baseline.is_finite() && optimized.is_finite());
                assert!(
                    (optimized - baseline).abs() <= 8.0 * f32::EPSILON * baseline.abs(),
                    "a={a}, powf={baseline}, sqrt={optimized}"
                );
            }
        }
    }

    #[test]
    fn eighth_root_threshold_matches_powf() {
        for eps in [f64::EPSILON / 4.0, 3.0 * f64::EPSILON] {
            let baseline = eps.powf(-1.0 / 8.0);
            let optimized = eps.sqrt().sqrt().sqrt().recip();
            assert!((optimized - baseline).abs() <= 2.0 * f64::EPSILON * baseline);
        }
        for eps in [f32::EPSILON / 4.0, 3.0 * f32::EPSILON] {
            let baseline = eps.powf(-1.0 / 8.0);
            let optimized = eps.sqrt().sqrt().sqrt().recip();
            assert!((optimized - baseline).abs() <= 2.0 * f32::EPSILON * baseline);
        }
    }
}
