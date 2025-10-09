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

use crate::{
    crate_util::{check, declare},
    polyeval, StrErr,
};

/// Computes [complete elliptic integral of the first kind](https://dlmf.nist.gov/19.2.E8).
/// ```text
///           π/2
///          ⌠          dθ
/// K(m)  =  │  _________________
///          │     _____________
///          ⌡   \╱ 1 - m sin²θ
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
/// ![Complete Elliptic Integral of the First Kind](https://github.com/p-sira/ellip/blob/main/figures/ellipk_plot.svg?raw=true)
///
/// [Interactive Plot](https://github.com/p-sira/ellip/blob/main/figures/ellipk_plot.html)
///
/// ## Special Cases
/// - K(0) = π/2
/// - K(1) = ∞
/// - K(-∞) = 0
///
/// # Related Functions
/// - [ellipk](crate::ellipk)(m) = [elliprf](crate::elliprf)(0, 1 - m, 1)
/// - [ellipf](crate::ellipf)(π/2, m) = [ellipk](crate::ellipk)(m)
///
/// # Examples
/// ```
/// use ellip::{ellipk, util::assert_close};
///
/// assert_close(ellipk(0.5).unwrap(), 1.8540746773013719, 1e-15);
/// ```
///
/// # References
/// - Maddock, John, Paul Bristow, Hubert Holin, and Xiaogang Zhang. “Boost Math Library: Special Functions - Elliptic Integrals.” Accessed April 17, 2025. <https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/ellint.html>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn ellipk<T: Float>(m: T) -> Result<T, StrErr> {
    match (m * 20.0).to_i8() {
        Some(0) | Some(1) => {
            let coeffs = [
                1.591003453790792180,
                0.416000743991786912,
                0.245791514264103415,
                0.179481482914906162,
                0.144556057087555150,
                0.123200993312427711,
                0.108938811574293531,
                0.098853409871592910,
                0.091439629201749751,
                0.085842591595413900,
                0.081541118718303215,
                0.078199656811256481910,
            ];
            Ok(polyeval(m - 0.05, &coeffs))
        }
        Some(2) | Some(3) => {
            let coeffs = [
                1.635256732264579992,
                0.471190626148732291,
                0.309728410831499587,
                0.252208311773135699,
                0.226725623219684650,
                0.215774446729585976,
                0.213108771877348910,
                0.216029124605188282,
                0.223255831633057896,
                0.234180501294209925,
                0.248557682972264071,
                0.266363809892617521,
            ];
            Ok(polyeval(m - 0.15, &coeffs))
        }
        Some(4) | Some(5) => {
            let coeffs = [
                1.685750354812596043,
                0.541731848613280329,
                0.401524438390690257,
                0.369642473420889090,
                0.376060715354583645,
                0.405235887085125919,
                0.453294381753999079,
                0.520518947651184205,
                0.609426039204995055,
                0.724263522282908870,
                0.871013847709812357,
                1.057652872753547036,
            ];
            Ok(polyeval(m - 0.25, &coeffs))
        }
        Some(6) | Some(7) => {
            let coeffs = [
                1.744350597225613243,
                0.634864275371935304,
                0.539842564164445538,
                0.571892705193787391,
                0.670295136265406100,
                0.832586590010977199,
                1.073857448247933265,
                1.422091460675497751,
                1.920387183402304829,
                2.632552548331654201,
                3.652109747319039160,
                5.115867135558865806,
                7.224080007363877411,
            ];
            Ok(polyeval(m - 0.35, &coeffs))
        }
        Some(8) | Some(9) => {
            let coeffs = [
                1.813883936816982644,
                0.763163245700557246,
                0.761928605321595831,
                0.951074653668427927,
                1.315180671703161215,
                1.928560693477410941,
                2.937509342531378755,
                4.594894405442878062,
                7.330071221881720772,
                11.87151259742530180,
                19.45851374822937738,
                32.20638657246426863,
                53.73749198700554656,
                90.27388602940998849,
            ];
            Ok(polyeval(m - 0.45, &coeffs))
        }
        Some(10) | Some(11) => {
            let coeffs = [
                1.898924910271553526,
                0.950521794618244435,
                1.151077589959015808,
                1.750239106986300540,
                2.952676812636875180,
                5.285800396121450889,
                9.832485716659979747,
                18.78714868327559562,
                36.61468615273698145,
                72.45292395127771801,
                145.1079577347069102,
                293.4786396308497026,
                598.3851815055010179,
                1228.420013075863451,
                2536.529755382764488,
            ];
            Ok(polyeval(m - 0.55, &coeffs))
        }
        Some(12) | Some(13) => {
            let coeffs = [
                2.007598398424376302,
                1.248457231212347337,
                1.926234657076479729,
                3.751289640087587680,
                8.119944554932045802,
                18.66572130873555361,
                44.60392484291437063,
                109.5092054309498377,
                274.2779548232413480,
                697.5598008606326163,
                1795.716014500247129,
                4668.381716790389910,
                12235.76246813664335,
                32290.17809718320818,
                85713.07608195964685,
                228672.1890493117096,
                612757.2711915852774,
            ];
            Ok(polyeval(m - 0.65, &coeffs))
        }
        Some(14) | Some(15) => {
            let coeffs = [
                2.156515647499643235,
                1.791805641849463243,
                3.826751287465713147,
                10.38672468363797208,
                31.40331405468070290,
                100.9237039498695416,
                337.3268282632272897,
                1158.707930567827917,
                4060.990742193632092,
                14454.00184034344795,
                52076.66107599404803,
                189493.6591462156887,
                695184.5762413896145,
                2567994.048255284686,
                9541921.966748386322,
                35634927.44218076174,
                133669298.4612040871,
                503352186.6866284541,
                1901975729.538660119,
                7208915015.330103756,
            ];
            Ok(polyeval(m - 0.75, &coeffs))
        }
        Some(16) => {
            let coeffs = [
                2.318122621712510589,
                2.616920150291232841,
                7.897935075731355823,
                30.50239715446672327,
                131.4869365523528456,
                602.9847637356491617,
                2877.024617809972641,
                14110.51991915180325,
                70621.44088156540229,
                358977.2665825309926,
                1847238.263723971684,
                9600515.416049214109,
                50307677.08502366879,
                265444188.6527127967,
                1408862325.028702687,
                7515687935.373774627,
            ];
            Ok(polyeval(m - 0.825, &coeffs))
        }
        Some(17) => {
            let coeffs = [
                2.473596173751343912,
                3.727624244118099310,
                15.60739303554930496,
                84.12850842805887747,
                506.9818197040613935,
                3252.277058145123644,
                21713.24241957434256,
                149037.0451890932766,
                1043999.331089990839,
                7427974.817042038995,
                53503839.67558661151,
                389249886.9948708474,
                2855288351.100810619,
                21090077038.76684053,
                156699833947.7902014,
                1170222242422.439893,
                8777948323668.937971,
                66101242752484.95041,
                499488053713388.7989,
                37859743397240299.20,
            ];
            Ok(polyeval(m - 0.875, &coeffs))
        }
        Some(_) => ellipk_precise(m),
        None => {
            check!(@nan, ellipk, [m]);
            if m == neg_inf!() {
                return Ok(0.0);
            }
            #[cfg(not(feature = "test_force_fail"))]
            if m > 1.0 {
                // Also handles inf
                return Err("ellipk: m must not be greater than 1.");
            }
            Err("ellipk: Unexpected error.")
        }
    }
}

#[inline]
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub(crate) fn ellipk_precise<T: Float>(m: T) -> Result<T, StrErr> {
    // Special cases: https://dlmf.nist.gov/19.6.E1
    if m >= 1.0 {
        if m == 1.0 {
            return Ok(inf!());
        }
        return Err("ellipk: m must not be greater than 1.");
    }

    Ok(ellipk_precise_unchecked(m))
}

/// Based on elliprf(1, 1-m, 0)
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub fn ellipk_precise_unchecked<T: Float>(m: T) -> T {
    declare!(mut [xn = T::one(), yn = (T::one() - m).sqrt(), t]);

    for _ in 0..MAX_ITERATION {
        if (xn - yn).abs() >= 2.7 * epsilon!() * xn.abs() {
            t = (xn * yn).sqrt();
            xn = (xn + yn) / 2.0;
            yn = t;
            continue;
        }
        break;
    }

    pi!() / (xn + yn)
}

const MAX_ITERATION: usize = 10;

#[cfg(not(feature = "test_force_fail"))]
#[cfg(all(test, not(feature = "no_std")))]
mod tests {
    use super::*;
    use crate::compare_test_data_boost;
    use crate::compare_test_data_wolfram;

    #[test]
    fn test_ellipk_boost() {
        compare_test_data_boost!("ellipk_data.txt", ellipk, 1, f64::EPSILON);
    }

    #[test]
    fn test_ellipk_wolfram() {
        compare_test_data_wolfram!("./tests/data/coverage", "ellipk_cov.csv", ellipk, 1, 2e-12);
        compare_test_data_wolfram!("ellipk_data.csv", ellipk, 1, 5e-15);
    }

    #[test]
    fn test_ellipk_special_cases() {
        use std::f64::{consts::FRAC_PI_2, INFINITY, NAN, NEG_INFINITY};
        // m = 0: K(0) = pi/2
        assert_eq!(ellipk(0.0).unwrap(), FRAC_PI_2);
        // m = 1: K(1) = inf
        assert_eq!(ellipk(1.0).unwrap(), INFINITY);
        // m < 0: should be valid, compare with reference value
        assert!(ellipk(-1.0).unwrap().is_finite());
        // m > 1: should return Err
        assert_eq!(ellipk(1.1), Err("ellipk: m must not be greater than 1."));
        // m = NaN: should return Err
        assert_eq!(ellipk(NAN), Err("ellipk: Arguments cannot be NAN."));
        // m = inf: should return Err
        assert_eq!(
            ellipk(INFINITY),
            Err("ellipk: m must not be greater than 1.")
        );
        // m = -inf: K(-inf) = 0
        assert_eq!(ellipk(NEG_INFINITY).unwrap(), 0.0);
    }
}

#[cfg(feature = "test_force_fail")]
crate::test_force_unreachable! {
    assert_eq!(ellipk(f64::INFINITY), Err("ellipk: Unexpected error."));
}
