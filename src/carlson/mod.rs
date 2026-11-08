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

    // Bit-exact reference values computed with the original powf-based
    // implementations at commit af8fe7d. Check the public RD/RG functions
    // as well as RJ(p=z), which dispatches to RD internally.
    #[test]
    fn paired_carlson_results_f64() {
        let cases = [
            (1.0, 0.5, 0.25, 0x40101723147c67fa, 0x3fe815e3e72b48d8),
            (
                1.1111111111111112,
                0.888888888888889,
                1.0,
                0x3ff00ae67aea8861,
                0x3feffc9f102d99c0,
            ),
            (0.0, 16.0, 1.0, 0x3fe6f811b0ed1ba2, 0x40012826e61f2f8c),
            (
                0.0,
                f64::from_bits(16.0_f64.to_bits() + 1),
                1.0,
                0x3fe6f811b0ed1ba4,
                0x40012826e61f2f8c,
            ),
            (
                1e-200,
                2e-200,
                3e-200,
                0x7e1bc218540fe564,
                0x2b339fa5c0e2fc50,
            ),
            (1e200, 2e200, 3e200, 0x0188e6012a86590c, 0x54b9a2fe925e9708),
            (
                167.001708984375,
                193.7735595703125,
                167.001708984375,
                0x3f3d05235dd745f3,
                0x402a8549db71c032,
            ),
            (2.0, 3.0, 0.25, 0x3ff9420c5d23267c, 0x3ff49b9e81137c08),
        ];
        for (x, y, z, rd_bits, rg_bits) in cases {
            let rd = f64::from_bits(rd_bits);
            let rg = f64::from_bits(rg_bits);
            for actual in [
                super::elliprd(x, y, z).unwrap(),
                super::elliprj(x, y, z, z).unwrap(),
            ] {
                assert!(actual.is_finite());
                assert!(
                    (actual - rd).abs() <= 8.0 * f64::EPSILON * rd.abs(),
                    "RD/RJ at ({x}, {y}, {z}): baseline={rd}, current={actual}"
                );
            }
            let actual = super::elliprg(x, y, z).unwrap();
            assert!(actual.is_finite());
            assert!(
                (actual - rg).abs() <= 8.0 * f64::EPSILON * rg.abs(),
                "RG at ({x}, {y}, {z}): baseline={rg}, current={actual}"
            );
        }
    }

    #[test]
    fn paired_carlson_results_f32() {
        let cases = [
            (1.0_f32, 0.5, 0.25, 0x4080b918, 0x3f40af1e),
            (1.3333334, 1.7777778, 1.6, 0x3f0134b5, 0x3fa04aee),
            (0.0, 16.0, 1.0, 0x3f37c08e, 0x40094137),
            (
                0.0,
                f32::from_bits(16.0_f32.to_bits() + 1),
                1.0,
                0x3f37c08e,
                0x40094137,
            ),
            (1e-20, 2e-20, 3e-20, 0x706aa1c7, 0x2f1a227c),
            (1e20, 2e20, 3e20, 0x0cbc84fd, 0x5050e438),
            (2.0, 3.0, 0.25, 0x3fca1064, 0x3fa4dcf4),
        ];
        for (x, y, z, rd_bits, rg_bits) in cases {
            let rd = f32::from_bits(rd_bits);
            let rg = f32::from_bits(rg_bits);
            for actual in [
                super::elliprd(x, y, z).unwrap(),
                super::elliprj(x, y, z, z).unwrap(),
            ] {
                assert!(actual.is_finite());
                assert!(
                    (actual - rd).abs() <= 8.0 * f32::EPSILON * rd.abs(),
                    "RD/RJ at ({x}, {y}, {z}): baseline={rd}, current={actual}"
                );
            }
            let actual = super::elliprg(x, y, z).unwrap();
            assert!(actual.is_finite());
            assert!(
                (actual - rg).abs() <= 8.0 * f32::EPSILON * rg.abs(),
                "RG at ({x}, {y}, {z}): baseline={rg}, current={actual}"
            );
        }
    }
}
