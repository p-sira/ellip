/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use std::f64::consts::PI;

// Reference: Derby and Olbert, “Cylindrical Magnets and Ideal Solenoids.”
/// Compute complete elliptic integral in Bulirsch form.
/// ```text
///                       π/2                                                   
///                      ⌠                   2              2
///                      ⎮            a ⋅ cos  (ϑ) + b ⋅ sin  (ϑ)               
/// cel(kc, p, a, b)  =  ⎮ ─────────────────────────────────────────────────────── ⋅ dϑ
///                      ⎮                              _________________________
///                      ⎮    2              2         ╱   2         2      2       
///                      ⌡ cos  (ϑ) + p ⋅ sin  (ϑ) ⋅ ╲╱ cos  (ϑ) + kc  ⋅ sin  (ϑ)    
///                     0                                                   
/// where kc ≠ 0, p ≠ 0
/// ```
pub fn cel(kc: f64, p: f64, a: f64, b: f64) -> Result<f64, &'static str> {
    if kc == 0.0 {
        return Err("cel: kc cannot be zero.");
    }

    if p == 0.0 {
        return Err("cel: p cannot be zero.");
    }

    let mut k = kc.abs();
    let mut aa;
    let mut pp;
    let mut bb;

    if p > 0.0 {
        aa = a;
        pp = p.sqrt();
        bb = b / pp;
    } else {
        let f = kc * kc;
        let q = (1.0 - f) * (b - a * p);
        let g = 1.0 - p;
        let h = f - p;
        pp = (h / g).sqrt();
        aa = (a - b) / g;
        bb = -q / (g * g * pp) + aa * pp;
    }

    let mut em = 1.0;
    let mut f = aa;
    aa += bb / pp;
    let mut g = k / pp;
    bb = 2.0 * (bb + f * g);
    pp += g;
    g = em;
    em += k;
    let mut kk = k;

    while (g - k).abs() > (g * ERRTOL) {
        k = 2.0 * kk.sqrt();
        kk = k * em;
        f = aa;
        aa += bb / pp;
        g = kk / pp;
        bb = 2.0 * (bb + f * g);
        pp += g;
        g = em;
        em += k;
    }

    Ok((PI / 2.0) * (bb + aa * em) / (em * (em + pp)))
}

const ERRTOL: f64 = 1e-6;

#[cfg(test)]
mod test {
    use super::*;
    use crate::{assert_close, ellipe, ellipk, test_util::linspace};

    /// Test using relationship with Legendre form.
    /// Reference: https://dlmf.nist.gov/19.2#iii
    #[test]
    fn test_cel() {
        fn test_kc(kc: f64) {
            let k = (1.0 - kc * kc).sqrt();
            let m = k * k;
            let ellipk = ellipk(m).unwrap();
            let ellipe = ellipe(m).unwrap();
            assert_close!(ellipk, cel(kc, 1.0, 1.0, 1.0).unwrap(), 1e-10);
            assert_close!(ellipe, cel(kc, 1.0, 1.0, kc * kc).unwrap(), 1e-10);
            assert_close!(
                (ellipe - kc * kc * ellipk) / m,
                cel(kc, 1.0, 1.0, 0.0).unwrap(),
                1e-10
            );
        }

        let linsp_neg = linspace(-1.0, -1e-3, 100);
        linsp_neg.iter().for_each(|kc| test_kc(*kc));
        let linsp_pos = linspace(1e-3, 1.0, 100);
        linsp_pos.iter().for_each(|kc| test_kc(*kc));
    }
}
