mod el1 {
    use ellip::*;

    #[test]
    fn test_el1() {
        compare_test_data_wolfram!("el1_data.csv", el1, 2, 5.0 * f64::EPSILON);
    }
}

mod el2 {
    use ellip::*;

    #[test]
    fn test_el2() {
        compare_test_data_wolfram!("el2_data.csv", el2, 4, 100.0 * f64::EPSILON);
    }
}

mod el3 {
    use ellip::*;

    #[test]
    fn test_el3() {
        compare_test_data_wolfram!("el3_data.csv", el3, 3, 3e-12);
        compare_test_data_wolfram!("el3_pv.csv", el3, 3, 5e-15);
    }
}

mod elliprc {
    use ellip::*;

    fn _elliprc(inp: &[f64]) -> f64 {
        elliprc(inp[0], inp[1]).unwrap()
    }

    #[test]
    fn test_elliprc() {
        compare_test_data_boost!("elliprc_data.txt", _elliprc, f64::EPSILON);
    }
}

mod elliprd {
    use ellip::*;

    fn _elliprd(inp: &[f64]) -> f64 {
        elliprd(inp[0], inp[1], inp[2]).unwrap()
    }

    #[test]
    fn test_elliprd() {
        compare_test_data_boost!("elliprd_data.txt", _elliprd, 4.8e-16);
    }

    #[test]
    fn test_elliprd_0xy() {
        compare_test_data_boost!("elliprd_0xy.txt", _elliprd, 5.9e-16);
    }

    #[test]
    fn test_elliprd_0yy() {
        compare_test_data_boost!("elliprd_0yy.txt", _elliprd, 2.6e-16);
    }

    #[test]
    fn test_elliprd_xxx() {
        compare_test_data_boost!("elliprd_xxx.txt", _elliprd, 2.3e-16);
    }

    #[test]
    fn test_elliprd_xxz() {
        compare_test_data_boost!("elliprd_xxz.txt", _elliprd, 7.9e-16);
    }

    #[test]
    fn test_elliprd_xyy() {
        compare_test_data_boost!("elliprd_xyy.txt", _elliprd, 3.7e-15);
    }
}

mod elliprf {
    use ellip::*;
    use itertools::Itertools;

    fn _elliprf(inp: &[f64]) -> f64 {
        let res = elliprf(inp[0], inp[1], inp[2]).unwrap();
        inp.iter().permutations(inp.len()).skip(1).for_each(|perm| {
            let actual = elliprf(*perm[0], *perm[1], *perm[2]).unwrap();
            assert_close!(res, actual, 6.5e-16);
        });
        res
    }

    #[test]
    fn test_elliprf() {
        compare_test_data_boost!("elliprf_data.txt", _elliprf, 4.5e-16);
    }

    #[test]
    fn test_elliprf_xxx() {
        compare_test_data_boost!("elliprf_xxx.txt", _elliprf, 2.3e-16);
    }

    #[test]
    fn test_elliprf_xy0() {
        compare_test_data_boost!("elliprf_xy0.txt", _elliprf, 4.2e-16);
    }

    #[test]
    fn test_elliprf_xyy() {
        compare_test_data_boost!("elliprf_xyy.txt", _elliprf, 5.3e-16);
    }

    #[test]
    fn test_elliprf_0yy() {
        compare_test_data_boost!("elliprf_0yy.txt", _elliprf, f64::EPSILON);
    }
}

mod elliprg {
    use ellip::*;
    use itertools::Itertools;

    fn _elliprg(inp: &[f64]) -> f64 {
        let res = elliprg(inp[0], inp[1], inp[2]).unwrap();
        inp.iter().permutations(inp.len()).skip(1).for_each(|perm| {
            let actual = elliprg(*perm[0], *perm[1], *perm[2]).unwrap();
            assert_close!(res, actual, 6.5e-16);
        });
        res
    }

    #[test]
    fn test_elliprg() {
        compare_test_data_boost!("elliprg_data.txt", _elliprg, 8.1e-16);
    }

    #[test]
    fn test_elliprg_wolfram() {
        compare_test_data_wolfram!("elliprg_data.csv", elliprg, 3, 1e-15);
    }

    #[test]
    fn test_elliprg_xxx() {
        compare_test_data_boost!("elliprg_xxx.txt", _elliprg, 2.4e-16);
    }

    #[test]
    fn test_elliprg_xy0() {
        compare_test_data_boost!("elliprg_xy0.txt", _elliprg, 4.4e-16);
    }

    #[test]
    fn test_elliprg_xyy() {
        compare_test_data_boost!("elliprg_xyy.txt", _elliprg, 5.4e-16);
    }

    #[test]
    fn test_elliprg_00x() {
        compare_test_data_boost!("elliprg_00x.txt", _elliprg, f64::EPSILON);
    }
}

mod elliprj {
    use ellip::*;
    use itertools::Itertools;

    fn _elliprj(inp: &[f64]) -> f64 {
        let res = elliprj(inp[0], inp[1], inp[2], inp[3]).unwrap();
        let (p, sym_params) = inp.split_last().unwrap();
        sym_params
            .iter()
            .permutations(sym_params.len())
            .skip(1)
            .for_each(|perm| {
                let actual = elliprj(*perm[0], *perm[1], *perm[2], *p).unwrap();
                assert_close!(res, actual, 3e-14);
            });
        res
    }

    #[test]
    fn test_elliprj_pv_rounding_f32() {
        // References rounded from 100-digit evaluations at exact f32 inputs.
        for (input, expected) in [
            (
                [2.0_f32 + 2.0 * f32::EPSILON, 2.0, 3.0, -1.0],
                1.2370731e-8_f32,
            ),
            ([2.0 - 2.0 * f32::EPSILON, 2.0, 3.0, -1.0], -1.2370734e-8),
            (
                [
                    2.0 + 2.0 * f32::EPSILON,
                    2.0 - 2.0 * f32::EPSILON,
                    3.0,
                    -1.0,
                ],
                -1.8478185e-15,
            ),
            ([0.8000001, 0.80000025, 1.2000004, -0.4000001], 5.277546e-9),
            ([0.8000001, 1.2000003, 0.80000037, -0.4000001], 3.0562223e-8),
            ([1.2, 0.80000025, 0.80000037, -0.4000001], 5.0569344e-8),
        ] {
            for scale in [2.0_f32.powi(-20), 1.0, 2.0_f32.powi(20)] {
                for perm in input[..3].iter().permutations(3) {
                    let actual = elliprj(
                        *perm[0] * scale,
                        *perm[1] * scale,
                        *perm[2] * scale,
                        input[3] * scale,
                    )
                    .unwrap()
                        * scale
                        * scale.sqrt();
                    assert!(
                        (actual / expected - 1.0).abs() <= 2.0 * f32::EPSILON,
                        "input={input:?}, scale={scale}, actual={actual}, expected={expected}"
                    );
                    assert_eq!(actual.is_sign_negative(), expected.is_sign_negative());
                }
            }
        }
        for scale in [2.0_f32.powi(-70), 1.0, 2.0_f32.powi(70)] {
            let actual = elliprj(scale, 2.0 * scale, 3.0 * scale, -4.0 * scale).unwrap()
                * scale
                * scale.sqrt();
            ellip::util::assert_close(actual, -0.2378677, 2.0 * f32::EPSILON);
            assert_eq!(
                elliprj(2.0 * scale, 2.0 * scale, 3.0 * scale, -scale),
                Ok(0.0)
            );
        }
        // Check the complete mu = one epsilon grid without an absolute floor.
        for case in ellip_dev_utils::parser::read_wolfram_data::<f32>(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/data/f32/wolfram/elliprj_pv.csv"
        ))
        .unwrap()
        {
            let actual = elliprj(
                case.inputs[0],
                case.inputs[1],
                case.inputs[2],
                case.inputs[3],
            )
            .unwrap();
            assert_eq!(actual, case.expected, "input={:?}", case.inputs);
        }
    }

    #[test]
    fn test_elliprj() {
        compare_test_data_boost!("elliprj_data.txt", _elliprj, 2.7e-14, atol: 5e-25);
        compare_test_data_boost!("elliprj_e2.txt", _elliprj, 4.8e-14, atol: 5e-25);
        compare_test_data_boost!("elliprj_e3.txt", _elliprj, 3.1e-15, atol: 5e-25);
        compare_test_data_boost!("elliprj_e4.txt", _elliprj, 2.2e-16, atol: 5e-25);
        compare_test_data_boost!("elliprj_zp.txt", _elliprj, 3.5e-15, atol: 5e-25);
    }

    #[test]
    fn test_elliprj_wolfram() {
        compare_test_data_wolfram!("elliprj_data.csv", elliprj, 4, 3.1e-15);
        compare_test_data_wolfram!("elliprj_pv.csv", elliprj, 4, 5e-14, atol: f64::EPSILON);
    }
}

mod ellipd {
    use ellip::*;

    #[test]
    fn test_ellipd() {
        compare_test_data_boost!("ellipd_data.txt", ellipd, 1, 2.9e-16);
    }

    #[test]
    fn test_ellipd_wolfram() {
        compare_test_data_wolfram!("ellipd_data.csv", ellipd, 1, 8e-16);
        compare_test_data_wolfram!("ellipd_neg.csv", ellipd, 1, 8e-16);
    }
}

mod ellipdinc {
    use ellip::*;

    #[test]
    fn test_ellipdinc() {
        compare_test_data_boost!("ellipdinc_data.txt", ellipdinc, 2, 6.4e-16);
    }

    #[test]
    fn test_ellipdinc_wolfram() {
        compare_test_data_wolfram!("ellipdinc_data.csv", ellipdinc, 2, 2e-15);
        compare_test_data_wolfram!("ellipdinc_neg.csv", ellipdinc, 2, 1e-15);
    }
}

mod ellipe {
    use ellip::*;

    #[test]
    fn test_ellipe() {
        compare_test_data_boost!("ellipe_data.txt", ellipe, 1, f64::EPSILON);
    }

    #[test]
    fn test_ellipe_wolfram() {
        compare_test_data_wolfram!("./tests/data/coverage", "ellipe_cov.csv", ellipe, 1, 7e-16);
    }
}

mod ellipeinc {
    use ellip::*;

    #[test]
    fn test_ellipeinc() {
        compare_test_data_boost!("ellipeinc_data.txt", ellipeinc, 2, 5e-16);
    }
}

mod ellipf {
    use ellip::*;

    #[test]
    fn test_ellipf() {
        compare_test_data_boost!("ellipf_data.txt", ellipf, 2, 5.1e-16);
    }
}

mod ellipk {
    use ellip::*;

    #[test]
    fn test_ellipk_boost() {
        compare_test_data_boost!("ellipk_data.txt", ellipk, 1, 2.8e-16);
    }

    #[test]
    fn test_ellipk_wolfram() {
        compare_test_data_wolfram!("./tests/data/coverage", "ellipk_cov.csv", ellipk, 1, 2e-12);
        compare_test_data_wolfram!("ellipk_data.csv", ellipk, 1, 5e-15);
    }
}

mod ellipke {
    use ellip::*;

    #[test]
    fn test_ellipke_wolfram() {
        fn ellipke_k(m: f64) -> Result<f64, StrErr> {
            Ok(ellipke(m)?.0)
        }
        fn ellipke_e(m: f64) -> Result<f64, StrErr> {
            Ok(ellipke(m)?.1)
        }

        compare_test_data_wolfram!("ellipk_data.csv", ellipke_k, 1, 5e-15);
        compare_test_data_wolfram!("ellipe_data.csv", ellipke_e, 1, 5e-15);
        compare_test_data_wolfram!("ellipk_neg.csv", ellipke_k, 1, 5e-15);
        compare_test_data_wolfram!("ellipe_neg.csv", ellipke_e, 1, 5e-15);
    }
}

mod ellippi {
    use ellip::*;

    #[test]
    fn test_ellippi() {
        compare_test_data_boost!("ellippi2_data_f64.txt", ellippi, 2, 6.6e-16);
    }

    #[test]
    fn test_ellippi_wolfram() {
        compare_test_data_wolfram!("ellippi_data.csv", ellippi, 2, 5e-15);
    }

    #[test]
    fn test_ellippi_wolfram_neg() {
        // A16 removes the stale dedicated dataset test because these cases now live here:
        // https://github.com/p-sira/ellip/pull/130
        // Includes the diagonal/near-diagonal cases generated by ellippi.wls.
        compare_test_data_wolfram!("ellippi_neg.csv", ellippi, 2, 1.3e-14);
    }

    #[test]
    fn test_ellippi_wolfram_pv() {
        compare_test_data_wolfram!("ellippi_pv.csv", ellippi, 2, 1e-14);
    }
}

mod ellippiinc {
    use ellip::*;

    fn ellippiinc_k(inp: &[f64]) -> f64 {
        ellippiinc(inp[1], inp[0], inp[2] * inp[2]).unwrap()
    }

    #[test]
    fn test_ellippiinc() {
        compare_test_data_boost!("ellippiinc_data.txt", ellippiinc_k, 2.1e-15);
        compare_test_data_boost!("ellippi3_large_data.txt", ellippiinc_k, 6e-15);
    }

    #[test]
    fn test_ellippiinc_wolfram() {
        compare_test_data_wolfram!("ellippiinc_data.csv", ellippiinc, 3, 5e-14);
    }

    #[test]
    fn test_ellippiinc_wolfram_neg() {
        compare_test_data_wolfram!("ellippiinc_neg.csv", ellippiinc, 3, 5e-14);
    }

    #[test]
    fn test_ellippiinc_wolfram_pv() {
        compare_test_data_wolfram!("ellippiinc_pv.csv", ellippiinc, 3, 3e-13);
    }

    #[test]
    fn test_ellippiinc_bulirsch_wolfram() {
        compare_test_data_wolfram!("ellippiinc_data.csv", ellippiinc_bulirsch, 3, 5e-14);
    }

    #[test]
    fn test_ellippiinc_bulirsch_wolfram_neg() {
        compare_test_data_wolfram!("ellippiinc_neg.csv", ellippiinc_bulirsch, 3, 5e-14);
    }

    #[test]
    fn test_ellippiinc_bulirsch_wolfram_pv() {
        compare_test_data_wolfram!("ellippiinc_pv.csv", ellippiinc_bulirsch, 3, 3e-13);
    }
}

mod heuman_lambda {
    use ellip::*;

    #[test]
    fn test_heuman_lambda() {
        compare_test_data_wolfram!("heuman_lambda_data.csv", heuman_lambda, 2, 3e-15);
    }
}

mod jacobi_zeta {
    use ellip::*;

    #[test]
    fn test_jacobi_zeta() {
        compare_test_data_wolfram!("jacobi_zeta_data.csv", jacobi_zeta, 2, 3e-15);
    }

    #[test]
    fn test_jacobi_neg() {
        compare_test_data_wolfram!("jacobi_zeta_neg.csv", jacobi_zeta, 2, 3e-15);
    }
}
