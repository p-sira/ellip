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
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  History:
//  XZ wrote the original of this file as part of the Google
//  Summer of Code 2006.  JM modified it to fit into the
//  Boost.Math conceptual framework better, and to ensure
//  that the code continues to work no matter how many digits
//  type T has.

use num_traits::Float;

use crate::{elliprf, polyeval};

/// Compute [complete elliptic integral of the first kind](https://dlmf.nist.gov/19.2.E8).
/// ```text
///           π/2
///          ⌠          dθ
/// K(m)  =  │  _________________
///          │     _____________
///          ⌡   \╱ 1 - m sin²θ
///         0
/// where m ≤ 1
/// ```
///
/// Note that some mathematical references use the parameter k for the function,
/// where k² = m.
///
/// # Examples
/// ```
/// use ellip::{ellipk, util::assert_close};
///
/// assert_close(ellipk(0.5).unwrap(), 1.8540746773013719, 1e-15);
/// ```
pub fn ellipk<T: Float>(m: T) -> Result<T, &'static str> {
    // If T is f128
    if max_val!() > num!(f64::MAX) {
        return ellipk_precise(m);
    }

    match (m * num!(20.0)).to_i32().unwrap() {
        0 | 1 => {
            let coeffs = [
                num!(1.591003453790792180),
                num!(0.416000743991786912),
                num!(0.245791514264103415),
                num!(0.179481482914906162),
                num!(0.144556057087555150),
                num!(0.123200993312427711),
                num!(0.108938811574293531),
                num!(0.098853409871592910),
                num!(0.091439629201749751),
                num!(0.085842591595413900),
                num!(0.081541118718303215),
                num!(0.078199656811256481910),
            ];
            Ok(polyeval(m - num!(0.05), &coeffs))
        }
        2 | 3 => {
            let coeffs = [
                num!(1.635256732264579992),
                num!(0.471190626148732291),
                num!(0.309728410831499587),
                num!(0.252208311773135699),
                num!(0.226725623219684650),
                num!(0.215774446729585976),
                num!(0.213108771877348910),
                num!(0.216029124605188282),
                num!(0.223255831633057896),
                num!(0.234180501294209925),
                num!(0.248557682972264071),
                num!(0.266363809892617521),
            ];
            Ok(polyeval(m - num!(0.15), &coeffs))
        }
        4 | 5 => {
            let coeffs = [
                num!(1.685750354812596043),
                num!(0.541731848613280329),
                num!(0.401524438390690257),
                num!(0.369642473420889090),
                num!(0.376060715354583645),
                num!(0.405235887085125919),
                num!(0.453294381753999079),
                num!(0.520518947651184205),
                num!(0.609426039204995055),
                num!(0.724263522282908870),
                num!(0.871013847709812357),
                num!(1.057652872753547036),
            ];
            Ok(polyeval(m - num!(0.25), &coeffs))
        }
        6 | 7 => {
            let coeffs = [
                num!(1.744350597225613243),
                num!(0.634864275371935304),
                num!(0.539842564164445538),
                num!(0.571892705193787391),
                num!(0.670295136265406100),
                num!(0.832586590010977199),
                num!(1.073857448247933265),
                num!(1.422091460675497751),
                num!(1.920387183402304829),
                num!(2.632552548331654201),
                num!(3.652109747319039160),
                num!(5.115867135558865806),
                num!(7.224080007363877411),
            ];
            Ok(polyeval(m - num!(0.35), &coeffs))
        }
        8 | 9 => {
            let coeffs = [
                num!(1.813883936816982644),
                num!(0.763163245700557246),
                num!(0.761928605321595831),
                num!(0.951074653668427927),
                num!(1.315180671703161215),
                num!(1.928560693477410941),
                num!(2.937509342531378755),
                num!(4.594894405442878062),
                num!(7.330071221881720772),
                num!(11.87151259742530180),
                num!(19.45851374822937738),
                num!(32.20638657246426863),
                num!(53.73749198700554656),
                num!(90.27388602940998849),
            ];
            Ok(polyeval(m - num!(0.45), &coeffs))
        }
        10 | 11 => {
            let coeffs = [
                num!(1.898924910271553526),
                num!(0.950521794618244435),
                num!(1.151077589959015808),
                num!(1.750239106986300540),
                num!(2.952676812636875180),
                num!(5.285800396121450889),
                num!(9.832485716659979747),
                num!(18.78714868327559562),
                num!(36.61468615273698145),
                num!(72.45292395127771801),
                num!(145.1079577347069102),
                num!(293.4786396308497026),
                num!(598.3851815055010179),
                num!(1228.420013075863451),
                num!(2536.529755382764488),
            ];
            Ok(polyeval(m - num!(0.55), &coeffs))
        }
        12 | 13 => {
            let coeffs = [
                num!(2.007598398424376302),
                num!(1.248457231212347337),
                num!(1.926234657076479729),
                num!(3.751289640087587680),
                num!(8.119944554932045802),
                num!(18.66572130873555361),
                num!(44.60392484291437063),
                num!(109.5092054309498377),
                num!(274.2779548232413480),
                num!(697.5598008606326163),
                num!(1795.716014500247129),
                num!(4668.381716790389910),
                num!(12235.76246813664335),
                num!(32290.17809718320818),
                num!(85713.07608195964685),
                num!(228672.1890493117096),
                num!(612757.2711915852774),
            ];
            Ok(polyeval(m - num!(0.65), &coeffs))
        }
        14 | 15 => {
            let coeffs = [
                num!(2.156515647499643235),
                num!(1.791805641849463243),
                num!(3.826751287465713147),
                num!(10.38672468363797208),
                num!(31.40331405468070290),
                num!(100.9237039498695416),
                num!(337.3268282632272897),
                num!(1158.707930567827917),
                num!(4060.990742193632092),
                num!(14454.00184034344795),
                num!(52076.66107599404803),
                num!(189493.6591462156887),
                num!(695184.5762413896145),
                num!(2567994.048255284686),
                num!(9541921.966748386322),
                num!(35634927.44218076174),
                num!(133669298.4612040871),
                num!(503352186.6866284541),
                num!(1901975729.538660119),
                num!(7208915015.330103756),
            ];
            Ok(polyeval(m - num!(0.75), &coeffs))
        }
        16 => {
            let coeffs = [
                num!(2.318122621712510589),
                num!(2.616920150291232841),
                num!(7.897935075731355823),
                num!(30.50239715446672327),
                num!(131.4869365523528456),
                num!(602.9847637356491617),
                num!(2877.024617809972641),
                num!(14110.51991915180325),
                num!(70621.44088156540229),
                num!(358977.2665825309926),
                num!(1847238.263723971684),
                num!(9600515.416049214109),
                num!(50307677.08502366879),
                num!(265444188.6527127967),
                num!(1408862325.028702687),
                num!(7515687935.373774627),
            ];
            Ok(polyeval(m - num!(0.825), &coeffs))
        }
        17 => {
            let coeffs = [
                num!(2.473596173751343912),
                num!(3.727624244118099310),
                num!(15.60739303554930496),
                num!(84.12850842805887747),
                num!(506.9818197040613935),
                num!(3252.277058145123644),
                num!(21713.24241957434256),
                num!(149037.0451890932766),
                num!(1043999.331089990839),
                num!(7427974.817042038995),
                num!(53503839.67558661151),
                num!(389249886.9948708474),
                num!(2855288351.100810619),
                num!(21090077038.76684053),
                num!(156699833947.7902014),
                num!(1170222242422.439893),
                num!(8777948323668.937971),
                num!(66101242752484.95041),
                num!(499488053713388.7989),
                num!(37859743397240299.20),
            ];
            Ok(polyeval(m - num!(0.875), &coeffs))
        }
        _ => ellipk_precise(m),
    }
}

#[inline]
pub(crate) fn ellipk_precise<T: Float>(m: T) -> Result<T, &'static str> {
    if m > one!() {
        return Err("ellipk: m must be less than 1.");
    }

    // Special cases: https://dlmf.nist.gov/19.6.E1
    if m == one!() {
        return Ok(inf!());
    }

    elliprf(zero!(), one!() - m, one!())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compare_test_data;

    fn ellipk_k(k: &[f64]) -> f64 {
        ellipk(k[0] * k[0]).unwrap()
    }

    #[test]
    fn test_ellipk() {
        compare_test_data!("./tests/data/boost/ellipk_data.txt", ellipk_k, f64::EPSILON);
    }
}
