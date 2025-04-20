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

use crate::{elliprg, polyeval};

/// Compute [complete elliptic integral of the second kind](<https://dlmf.nist.gov/19>.2.E8).
///
/// ```text
///           π/2
///          ⌠     ___________
/// E(m)  =  │ \╱ 1 - m sin²θ  dθ
///          ⌡
///         0
/// where m ≤ 1
/// ```
///
/// Note that some mathematical references use the parameter k for the function,
/// where k² = m.
///
/// # Examples
/// ```
/// use ellip::{ellipe, util::assert_close};
///
/// assert_close(ellipe(0.5).unwrap(), 1.3506438810476755, 1e-15);
/// ```
pub fn ellipe<T: Float>(m: T) -> Result<T, &'static str> {
    // If T is f128
    if max_val!() > num!(f64::MAX) {
        return ellipe_precise(m);
    }

    // Note: this function allows both negative m and positive m less than 1.
    // Negative m: Abramowitz & Stegun, 1972
    let mut m = m;
    let mut c = one!();
    while m < zero!() {
        c = c * (one!() - m).sqrt();
        m = m / (m - one!());
    }

    Ok(c * _ellipe(m)?)
}

#[inline]
fn _ellipe<T: Float>(m: T) -> Result<T, &'static str> {
    match (m * num!(20.0)).to_i32().unwrap() {
        0 | 1 => {
            let coeffs = [
                num!(1.550973351780472328),
                -num!(0.400301020103198524),
                -num!(0.078498619442941939),
                -num!(0.034318853117591992),
                -num!(0.019718043317365499),
                -num!(0.013059507731993309),
                -num!(0.009442372874146547),
                -num!(0.007246728512402157),
                -num!(0.005807424012956090),
                -num!(0.004809187786009338),
                -num!(0.004086399233255150),
            ];
            Ok(polyeval(m - num!(0.05), &coeffs))
        }
        2 | 3 => {
            let coeffs = [
                num!(1.510121832092819728),
                -num!(0.417116333905867549),
                -num!(0.090123820404774569),
                -num!(0.043729944019084312),
                -num!(0.027965493064761785),
                -num!(0.020644781177568105),
                -num!(0.016650786739707238),
                -num!(0.014261960828842520),
                -num!(0.012759847429264803),
                -num!(0.011799303775587354),
                -num!(0.011197445703074968),
            ];
            Ok(polyeval(m - num!(0.15), &coeffs))
        }
        4 | 5 => {
            let coeffs = [
                num!(1.467462209339427155),
                -num!(0.436576290946337775),
                -num!(0.105155557666942554),
                -num!(0.057371843593241730),
                -num!(0.041391627727340220),
                -num!(0.034527728505280841),
                -num!(0.031495443512532783),
                -num!(0.030527000890325277),
                -num!(0.030916984019238900),
                -num!(0.032371395314758122),
                -num!(0.034789960386404158),
            ];
            Ok(polyeval(m - num!(0.25), &coeffs))
        }
        6 | 7 => {
            let coeffs = [
                num!(1.422691133490879171),
                -num!(0.459513519621048674),
                -num!(0.125250539822061878),
                -num!(0.078138545094409477),
                -num!(0.064714278472050002),
                -num!(0.062084339131730311),
                -num!(0.065197032815572477),
                -num!(0.072793895362578779),
                -num!(0.084959075171781003),
                -num!(0.102539850131045997),
                -num!(0.127053585157696036),
                -num!(0.160791120691274606),
            ];
            Ok(polyeval(m - num!(0.35), &coeffs))
        }
        8 | 9 => {
            let coeffs = [
                num!(1.375401971871116291),
                -num!(0.487202183273184837),
                -num!(0.153311701348540228),
                -num!(0.111849444917027833),
                -num!(0.108840952523135768),
                -num!(0.122954223120269076),
                -num!(0.152217163962035047),
                -num!(0.200495323642697339),
                -num!(0.276174333067751758),
                -num!(0.393513114304375851),
                -num!(0.575754406027879147),
                -num!(0.860523235727239756),
                -num!(1.308833205758540162),
            ];
            Ok(polyeval(m - num!(0.45), &coeffs))
        }
        10 | 11 => {
            let coeffs = [
                num!(1.325024497958230082),
                -num!(0.521727647557566767),
                -num!(0.194906430482126213),
                -num!(0.171623726822011264),
                -num!(0.202754652926419141),
                -num!(0.278798953118534762),
                -num!(0.420698457281005762),
                -num!(0.675948400853106021),
                -num!(1.136343121839229244),
                -num!(1.976721143954398261),
                -num!(3.531696773095722506),
                -num!(6.446753640156048150),
                -num!(11.97703130208884026),
            ];
            Ok(polyeval(m - num!(0.55), &coeffs))
        }
        12 | 13 => {
            let coeffs = [
                num!(1.270707479650149744),
                -num!(0.566839168287866583),
                -num!(0.262160793432492598),
                -num!(0.292244173533077419),
                -num!(0.440397840850423189),
                -num!(0.774947641381397458),
                -num!(1.498870837987561088),
                -num!(3.089708310445186667),
                -num!(6.667595903381001064),
                -num!(14.89436036517319078),
                -num!(34.18120574251449024),
                -num!(80.15895841905397306),
                -num!(191.3489480762984920),
                -num!(463.5938853480342030),
                -num!(1137.380822169360061),
            ];
            Ok(polyeval(m - num!(0.65), &coeffs))
        }
        14 | 15 => {
            let coeffs = [
                num!(1.211056027568459525),
                -num!(0.630306413287455807),
                -num!(0.387166409520669145),
                -num!(0.592278235311934603),
                -num!(1.237555584513049844),
                -num!(3.032056661745247199),
                -num!(8.181688221573590762),
                -num!(23.55507217389693250),
                -num!(71.04099935893064956),
                -num!(221.8796853192349888),
                -num!(712.1364793277635425),
                -num!(2336.125331440396407),
                -num!(7801.945954775964673),
                -num!(26448.19586059191933),
                -num!(90799.48341621365251),
                -num!(315126.0406449163424),
                -num!(1104011.344311591159),
            ];
            Ok(polyeval(m - num!(0.75), &coeffs))
        }
        16 => {
            let coeffs = [
                num!(1.161307152196282836),
                -num!(0.701100284555289548),
                -num!(0.580551474465437362),
                -num!(1.243693061077786614),
                -num!(3.679383613496634879),
                -num!(12.81590924337895775),
                -num!(49.25672530759985272),
                -num!(202.1818735434090269),
                -num!(869.8602699308701437),
                -num!(3877.005847313289571),
                -num!(17761.70710170939814),
                -num!(83182.69029154232061),
                -num!(396650.4505013548170),
                -num!(1920033.413682634405),
            ];
            Ok(polyeval(m - num!(0.825), &coeffs))
        }
        17 => {
            let coeffs = [
                num!(1.124617325119752213),
                -num!(0.770845056360909542),
                -num!(0.844794053644911362),
                -num!(2.490097309450394453),
                -num!(10.23971741154384360),
                -num!(49.74900546551479866),
                -num!(267.0986675195705196),
                -num!(1532.665883825229947),
                -num!(9222.313478526091951),
                -num!(57502.51612140314030),
                -num!(368596.1167416106063),
                -num!(2415611.088701091428),
                -num!(16120097.81581656797),
                -num!(109209938.5203089915),
                -num!(749380758.1942496220),
                -num!(5198725846.725541393),
                -num!(36409256888.12139973),
            ];
            Ok(polyeval(m - num!(0.875), &coeffs))
        }
        _ => ellipe_precise(m),
    }
}

#[inline]
fn ellipe_precise<T: Float>(m: T) -> Result<T, &'static str> {
    if m > one!() {
        return Err("ellipe: m must be less than 1.");
    }

    // Special cases: <https://dlmf.nist.gov/19>.6.E1
    if m == one!() {
        return Ok(one!());
    }

    Ok(two!() * elliprg(zero!(), one!() - m, one!())?)
}

#[cfg(test)]
mod tests {
    use core::f64;

    use super::*;
    use crate::compare_test_data;

    fn ellipe_k(k: &[f64]) -> f64 {
        ellipe(k[0] * k[0]).unwrap()
    }

    #[test]
    fn test_ellipe() {
        compare_test_data!("./tests/data/boost/ellipe_data.txt", ellipe_k, f64::EPSILON);
    }
}
