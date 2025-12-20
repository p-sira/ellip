/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Boost Math, see LICENSE in this directory.
 */

// Original header from Boost Math
//  Copyright (c) 2006 Xiaogang Zhang
//  Copyright (c) 2006 John Maddock
//  Copyright (c) 2024 Matt Borland
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0.

use num_traits::Float;

use crate::{carlson::elliprg_unchecked, crate_util::check, polyeval, StrErr};

/// Computes [complete elliptic integral of the second kind](https://dlmf.nist.gov/19.2.E8).
/// ```text
///           π/2
///          ⌠     ___________
/// E(m)  =  │ \╱ 1 - m sin²θ  dθ
///          ⌡
///         0
/// ```
///
/// ## Parameters
/// - m: elliptic parameter. m ∈ ℝ, m ≤ 1.
///
/// The elliptic modulus (k) is also frequently used instead of the parameter (m), where k² = m.
///
/// ## Domain
/// - Returns error if m > 1.
///
/// ## Graph
/// ![Complete Elliptic Integral of the Second Kind](https://github.com/p-sira/ellip/blob/main/figures/ellipe.svg?raw=true)
///
/// [Interactive Plot](https://p-sira.github.io/ellippy/_static/figures/ellipe.html)
///
/// ## Special Cases
/// - E(0) = π/2
/// - E(1) = 1
/// - E(-∞) = ∞
///
/// # Related Functions
/// - [ellipe](crate::ellipe)(m) = 2 [elliprg](crate::elliprg)(0, 1 - m, 1)
/// - [ellipeinc](crate::ellipeinc)(π/2, m) = [ellipe](crate::ellipe)(m)
///
/// # Examples
/// ```
/// use ellip::{ellipe, util::assert_close};
///
/// assert_close(ellipe(0.5).unwrap(), 1.3506438810476755, 1e-15);
/// ```
///
/// # References
/// - Maddock, John, Paul Bristow, Hubert Holin, and Xiaogang Zhang. “Boost Math Library: Special Functions - Elliptic Integrals.” Accessed April 17, 2025. <https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/ellint.html>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
/// - Abramowitz, Milton, and Irene A. Stegun. Handbook of Mathematical Functions: With Formulas, Graphs and Mathematical Tables. Unabridged, Unaltered and corr. Republ. of the 1964 ed. With Conference on mathematical tables, National science foundation, and Massachusetts institute of technology. Dover Books on Advanced Mathematics. Dover publ, 1972.
/// - The SciPy community. “Scipy.Special.Ellipe — SciPy v1.16.0 Manual.” Accessed July 28, 2025. <https://docs.scipy.org/doc/scipy-1.16.0/reference/generated/scipy.special.ellipe.html>.
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn ellipe<T: Float>(m: T) -> Result<T, StrErr> {
    let mut m = m;
    let mut c = 1.0;
    if m < 0.0 {
        if m == neg_inf!() {
            return Ok(inf!());
        }
        // Negative m: Abramowitz & Stegun, 1972
        c = c * (1.0 - m).sqrt();
        m = m / (m - 1.0);
    }

    Ok(c * _ellipe(m)?)
}

#[inline]
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
fn _ellipe<T: Float>(m: T) -> Result<T, StrErr> {
    match (m * 20.0).to_i64() {
        Some(0) | Some(1) => {
            let coeffs = [
                1.550973351780472328,
                -0.400301020103198524,
                -0.078498619442941939,
                -0.034318853117591992,
                -0.019718043317365499,
                -0.013059507731993309,
                -0.009442372874146547,
                -0.007246728512402157,
                -0.005807424012956090,
                -0.004809187786009338,
                -0.004086399233255150,
            ];
            Ok(polyeval(m - 0.05, &coeffs))
        }
        Some(2) | Some(3) => {
            let coeffs = [
                1.510121832092819728,
                -0.417116333905867549,
                -0.090123820404774569,
                -0.043729944019084312,
                -0.027965493064761785,
                -0.020644781177568105,
                -0.016650786739707238,
                -0.014261960828842520,
                -0.012759847429264803,
                -0.011799303775587354,
                -0.011197445703074968,
            ];
            Ok(polyeval(m - 0.15, &coeffs))
        }
        Some(4) | Some(5) => {
            let coeffs = [
                1.467462209339427155,
                -0.436576290946337775,
                -0.105155557666942554,
                -0.057371843593241730,
                -0.041391627727340220,
                -0.034527728505280841,
                -0.031495443512532783,
                -0.030527000890325277,
                -0.030916984019238900,
                -0.032371395314758122,
                -0.034789960386404158,
            ];
            Ok(polyeval(m - 0.25, &coeffs))
        }
        Some(6) | Some(7) => {
            let coeffs = [
                1.422691133490879171,
                -0.459513519621048674,
                -0.125250539822061878,
                -0.078138545094409477,
                -0.064714278472050002,
                -0.062084339131730311,
                -0.065197032815572477,
                -0.072793895362578779,
                -0.084959075171781003,
                -0.102539850131045997,
                -0.127053585157696036,
                -0.160791120691274606,
            ];
            Ok(polyeval(m - 0.35, &coeffs))
        }
        Some(8) | Some(9) => {
            let coeffs = [
                1.375401971871116291,
                -0.487202183273184837,
                -0.153311701348540228,
                -0.111849444917027833,
                -0.108840952523135768,
                -0.122954223120269076,
                -0.152217163962035047,
                -0.200495323642697339,
                -0.276174333067751758,
                -0.393513114304375851,
                -0.575754406027879147,
                -0.860523235727239756,
                -1.308833205758540162,
            ];
            Ok(polyeval(m - 0.45, &coeffs))
        }
        Some(10) | Some(11) => {
            let coeffs = [
                1.325024497958230082,
                -0.521727647557566767,
                -0.194906430482126213,
                -0.171623726822011264,
                -0.202754652926419141,
                -0.278798953118534762,
                -0.420698457281005762,
                -0.675948400853106021,
                -1.136343121839229244,
                -1.976721143954398261,
                -3.531696773095722506,
                -6.446753640156048150,
                -11.97703130208884026,
            ];
            Ok(polyeval(m - 0.55, &coeffs))
        }
        Some(12) | Some(13) => {
            let coeffs = [
                1.270707479650149744,
                -0.566839168287866583,
                -0.262160793432492598,
                -0.292244173533077419,
                -0.440397840850423189,
                -0.774947641381397458,
                -1.498870837987561088,
                -3.089708310445186667,
                -6.667595903381001064,
                -14.89436036517319078,
                -34.18120574251449024,
                -80.15895841905397306,
                -191.3489480762984920,
                -463.5938853480342030,
                -1137.380822169360061,
            ];
            Ok(polyeval(m - 0.65, &coeffs))
        }
        Some(14) | Some(15) => {
            let coeffs = [
                1.211056027568459525,
                -0.630306413287455807,
                -0.387166409520669145,
                -0.592278235311934603,
                -1.237555584513049844,
                -3.032056661745247199,
                -8.181688221573590762,
                -23.55507217389693250,
                -71.04099935893064956,
                -221.8796853192349888,
                -712.1364793277635425,
                -2336.125331440396407,
                -7801.945954775964673,
                -26448.19586059191933,
                -90799.48341621365251,
                -315126.0406449163424,
                -1104011.344311591159,
            ];
            Ok(polyeval(m - 0.75, &coeffs))
        }
        Some(16) => {
            let coeffs = [
                1.161307152196282836,
                -0.701100284555289548,
                -0.580551474465437362,
                -1.243693061077786614,
                -3.679383613496634879,
                -12.81590924337895775,
                -49.25672530759985272,
                -202.1818735434090269,
                -869.8602699308701437,
                -3877.005847313289571,
                -17761.70710170939814,
                -83182.69029154232061,
                -396650.4505013548170,
                -1920033.413682634405,
            ];
            Ok(polyeval(m - 0.825, &coeffs))
        }
        Some(17) => {
            let coeffs = [
                1.124617325119752213,
                -0.770845056360909542,
                -0.844794053644911362,
                -2.490097309450394453,
                -10.23971741154384360,
                -49.74900546551479866,
                -267.0986675195705196,
                -1532.665883825229947,
                -9222.313478526091951,
                -57502.51612140314030,
                -368596.1167416106063,
                -2415611.088701091428,
                -16120097.81581656797,
                -109209938.5203089915,
                -749380758.1942496220,
                -5198725846.725541393,
                -36409256888.12139973,
            ];
            Ok(polyeval(m - 0.875, &coeffs))
        }
        Some(_) => ellipe_precise(m),
        None => {
            check!(@nan, ellipe, [m]);
            #[cfg(not(feature = "test_force_fail"))]
            if m > 1.0 {
                // Infinity cases
                return Err("ellipe: m must not be greater than 1.");
            }
            Err("ellipe: Unexpected error.")
        }
    }
}

#[inline]
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
fn ellipe_precise<T: Float>(m: T) -> Result<T, StrErr> {
    // Special cases: https://dlmf.nist.gov/19.6.E1
    if m >= 1.0 {
        if m == 1.0 {
            return Ok(1.0);
        }
        return Err("ellipe: m must not be greater than 1.");
    }

    Ok(2.0 * elliprg_unchecked(0.0, 1.0 - m, 1.0))
}

#[cfg(not(feature = "test_force_fail"))]
#[cfg(all(test, not(feature = "no_std")))]
mod tests {
    use core::f64;

    use super::*;
    use crate::{compare_test_data_boost, compare_test_data_wolfram};

    #[test]
    fn test_ellipe() {
        compare_test_data_boost!("ellipe_data.txt", ellipe, 1, f64::EPSILON);
    }

    #[test]
    fn test_ellipe_wolfram() {
        compare_test_data_wolfram!("./tests/data/coverage", "ellipe_cov.csv", ellipe, 1, 7e-16);
    }

    #[test]
    fn test_ellipe_special_cases() {
        use std::f64::{consts::FRAC_PI_2, INFINITY, NAN, NEG_INFINITY};
        // m > 1: should return Err
        assert_eq!(ellipe(1.1), Err("ellipe: m must not be greater than 1."));
        // m = 0: E(0) = pi/2
        assert_eq!(ellipe(0.0).unwrap(), FRAC_PI_2);
        // m = 1: E(1) = 1
        assert_eq!(ellipe(1.0).unwrap(), 1.0);
        // m < 0: should be valid, compare with reference value
        assert!(ellipe(-1.0).unwrap().is_finite());
        // m = NaN: should return Err
        assert_eq!(ellipe(NAN), Err("ellipe: Arguments cannot be NAN."));
        // m = inf: should return Err
        assert_eq!(
            ellipe(INFINITY),
            Err("ellipe: m must not be greater than 1.")
        );
        // m = -inf: E(-inf) = inf
        assert_eq!(ellipe(NEG_INFINITY).unwrap(), INFINITY);
    }
}

#[cfg(feature = "test_force_fail")]
crate::test_force_unreachable! {
    assert_eq!(ellipe(f64::INFINITY), Err("ellipe: Unexpected error."));
}
