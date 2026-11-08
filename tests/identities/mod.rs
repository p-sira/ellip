mod cel {
    use crate::common::linspace;
    use ellip::*;
    use itertools::iproduct;

    /// Test using relationship with Legendre form.
    /// Reference: https://dlmf.nist.gov/19.2#iii
    #[test]
    fn test_cel() {
        fn _test(kc: f64, p: f64) {
            let m = 1.0 - kc * kc;
            let ellipk = ellipk(m).unwrap();
            let ellipe = ellipe(m).unwrap();

            // cel precision is low for K cases
            assert_close! {ellipk, cel(kc, 1.0, 1.0, 1.0).unwrap(), 2e-12};
            assert_close! {ellipe, cel(kc, 1.0, 1.0, kc * kc).unwrap(), 1e-15};
            assert_close! {(ellipe - kc * kc * ellipk) / m, cel(kc, 1.0, 1.0, 0.0).unwrap(), 8.5e-14};

            // Bulirsch, “Numerical Calculation of Elliptic Integrals and Elliptic Functions III”
            let n = 1.0 - p;
            let ellippi = ellippi(n, m).unwrap();

            assert_close! {(ellipk - ellipe) / m, cel(kc, 1.0, 0.0, 1.0).unwrap(), 6e-12};
            // cel precision is very low for PI cases
            assert_close! {ellippi, cel(kc, p, 1.0, 1.0).unwrap(), 3.5e-12};
            assert_close! {(ellippi - ellipk) / (1.0 - p), cel(kc, p, 0.0, 1.0).unwrap(), 3.5e-12};
        }

        let linsp_kc = [
            linspace(-1.0 + 1e-3, -1e-3, 100),
            linspace(1e-3, 1.0 - 1e-3, 100),
        ]
        .concat();
        let linsp_p = linspace(1e-3, 1.0 - 1e-3, 10);

        iproduct!(linsp_kc, linsp_p).for_each(|(kc, p)| _test(kc, p));

        // Data from Bulirsch, “Numerical Calculation of Elliptic Integrals and Elliptic Functions III”
        assert_close! {cel(1e-1, 4.1, 1.2, 1.1).unwrap(), 1.5464442694017956, 5e-16};
        assert_close! {cel(1e-1, -4.1, 1.2, 1.1).unwrap(), -6.7687378198360556e-1, 5e-16};
    }
}

mod cel1 {
    use crate::common::linspace;
    use ellip::*;

    #[test]
    fn test_cel1() {
        fn test_kc(kc: f64) {
            let m = 1.0 - kc * kc;
            assert_close!(ellipk(m).unwrap(), cel1(kc).unwrap(), 2e-12);
        }

        let linsp_neg = linspace(-1.0, -1e-3, 100);
        linsp_neg.iter().for_each(|kc| test_kc(*kc));
        let linsp_pos = linspace(1e-3, 1.0, 100);
        linsp_pos.iter().for_each(|kc| test_kc(*kc));
    }
}

mod cel2 {
    use crate::common::linspace;
    use ellip::*;

    #[test]
    fn test_cel2() {
        fn test_kc(kc: f64) {
            let m = 1.0 - kc * kc;
            assert_close!(ellipe(m).unwrap(), cel2(kc, 1.0, kc * kc).unwrap(), 7.7e-16);
        }

        let linsp_neg = linspace(-1.0, -1e-3, 100);
        linsp_neg.iter().for_each(|kc| test_kc(*kc));
        let linsp_pos = linspace(1e-3, 1.0, 100);
        linsp_pos.iter().for_each(|kc| test_kc(*kc));
    }
}

mod cel3 {
    use crate::common::linspace;
    use ellip::*;

    #[test]
    fn test_cel3() {
        for kc in linspace(0.05, 5.0, 20) {
            // p = 1 matches cel1
            assert_close! {cel3(kc, 1.0).unwrap(), cel1(kc).unwrap(), 1e-15};

            for p in linspace(0.05, 5.0, 20) {
                let actual = cel3(kc, p).unwrap();
                let expected = cel(kc, p, 1.0, 1.0).unwrap();
                assert_close! {actual, expected, 1e-15};
            }

            for p in linspace(-5.0, -0.05, 20) {
                let actual = cel3(kc, p).unwrap();
                let expected = cel(kc, p, 1.0, 1.0).unwrap();
                assert_close! {actual, expected, 1e-14};
            }
        }
    }

    #[test]
    fn test_cel3_independent_reference() {
        use ellip::{ellipk, ellippi};

        // For kc in (0, 1), kc^2 = 1 - m => m = 1 - kc^2.
        // With n = 1 - p: cel3(kc, p) = ellippi(1 - p, 1 - kc^2).
        // ellippi uses Carlson RF and RJ, an independent recurrence.
        for kc in linspace(0.1, 0.9, 10) {
            let m = 1.0 - kc * kc;
            // p = 1 matches ellipk(m)
            assert_close! {cel3(kc, 1.0).unwrap(), ellipk(m).unwrap(), 1e-14};

            // Positive p (n < 1 and n > 1):
            for p in [0.2, 0.5, 0.8, 1.5, 2.0] {
                let actual = cel3(kc, p).unwrap();
                let expected = ellippi(1.0 - p, m).unwrap();
                assert_close! {actual, expected, 1e-12};
            }
        }
    }
}

mod carlson {
    use ellip::{elliprd, elliprg, elliprj};

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
            for actual in [elliprd(x, y, z).unwrap(), elliprj(x, y, z, z).unwrap()] {
                assert!(actual.is_finite());
                assert!(
                    (actual - rd).abs() <= 8.0 * f64::EPSILON * rd.abs(),
                    "RD/RJ at ({x}, {y}, {z}): baseline={rd}, current={actual}"
                );
            }
            let actual = elliprg(x, y, z).unwrap();
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
            for actual in [elliprd(x, y, z).unwrap(), elliprj(x, y, z, z).unwrap()] {
                assert!(actual.is_finite());
                assert!(
                    (actual - rd).abs() <= 8.0 * f32::EPSILON * rd.abs(),
                    "RD/RJ at ({x}, {y}, {z}): baseline={rd}, current={actual}"
                );
            }
            let actual = elliprg(x, y, z).unwrap();
            assert!(actual.is_finite());
            assert!(
                (actual - rg).abs() <= 8.0 * f32::EPSILON * rg.abs(),
                "RG at ({x}, {y}, {z}): baseline={rg}, current={actual}"
            );
        }
    }
}

mod ellipke {
    use crate::common::linspace;
    use ellip::util::assert_close;
    use ellip::*;

    #[test]
    fn test_ellipke_matches_individual() {
        // Test positive range [0, 0.999]
        for m in linspace(0.0, 0.999, 100) {
            let (k, e) = ellipke(m).unwrap();
            let k_ref = ellipk(m).unwrap();
            let e_ref = ellipe(m).unwrap();
            assert_close(k, k_ref, 1e-14);
            assert_close(e, e_ref, 1e-14);
        }

        // Test negative range [-10.0, 0.0]
        for m in linspace(-10.0, -0.001, 100) {
            let (k, e) = ellipke(m).unwrap();
            let k_ref = ellipk(m).unwrap();
            let e_ref = ellipe(m).unwrap();
            assert_close(k, k_ref, 1e-14);
            assert_close(e, e_ref, 1e-14);
        }

        // Test extreme negative range
        for &m in &[-1e4, -1e8, -1e12, -1e16, -1e50] {
            let (k, e) = ellipke(m).unwrap();
            let k_ref = ellipk(m).unwrap();
            let e_ref = ellipe(m).unwrap();
            assert_close(k, k_ref, 1e-13);
            assert_close(e, e_ref, 1e-13);
        }
    }
}
