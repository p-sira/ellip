/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use std::f64::consts::FRAC_PI_2;

use num_traits::Float;

// Reference: Derby and Olbert, “Cylindrical Magnets and Ideal Solenoids.”
/// Compute [complete elliptic integral in Bulirsch form](https://dlmf.nist.gov/19.2#iii).
/// ```text
///                       π/2                                                   
///                      ⌠            a cos²(ϑ) + b sin²(ϑ)               
/// cel(kc, p, a, b)  =  ⎮ ─────────────────────────────────────────────── ⋅ dϑ
///                      ⎮                         ______________________
///                      ⌡ (cos²(ϑ) + p sin²(ϑ)) ╲╱ cos²(ϑ) + kc² sin²(ϑ)    
///                     0                                                   
/// where kc ≠ 0, p ≠ 0
/// ```
pub fn cel<T: Float>(kc: T, p: T, a: T, b: T) -> Result<T, &'static str> {
    if kc == T::zero() {
        return Err("cel: kc cannot be zero.");
    }

    if p == T::zero() {
        return Err("cel: p cannot be zero.");
    }

    Ok(_cel(kc, p, a, b))
}

/// Unchecked version of [cel].
pub fn _cel<T: Float>(kc: T, p: T, a: T, b: T) -> T {
    let mut k = kc.abs();
    let mut aa;
    let mut pp;
    let mut bb;

    if p > T::zero() {
        aa = a;
        pp = p.sqrt();
        bb = b / pp;
    } else {
        let f = kc * kc;
        let q: T = (T::one() - f) * (b - a * p);
        let g = T::one() - p;
        let h = f - p;
        pp = (h / g).sqrt();
        aa = (a - b) / g;
        bb = -q / (g * g * pp) + aa * pp;
    }

    let mut em = T::one();
    let mut f = aa;
    aa = aa + bb / pp;
    let mut g = k / pp;
    bb = T::from(2.0).unwrap() * (bb + f * g);
    pp = pp + g;
    g = em;
    em = em + k;
    let mut kk = k;

    while (g - k).abs() > (g * T::from(1e-6).unwrap()) {
        k = T::from(2.0).unwrap() * kk.sqrt();
        kk = k * em;
        f = aa;
        aa = aa + bb / pp;
        g = kk / pp;
        bb = T::from(2.0).unwrap() * (bb + f * g);
        pp = pp + g;
        g = em;
        em = em + k;
    }

    T::from(FRAC_PI_2).unwrap() * (bb + aa * em) / (em * (em + pp))
}

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
