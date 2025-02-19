/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use num_traits::Float;

use super::BulirschConst;

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
pub fn cel<T: Float + BulirschConst>(kc: T, p: T, a: T, b: T) -> Result<T, &'static str> {
    if kc == zero!() {
        return Err("cel: kc cannot be zero.");
    }

    if p == zero!() {
        return Err("cel: p cannot be zero.");
    }

    let mut k = kc.abs();
    let mut aa;
    let mut pp;
    let mut bb: T;

    if p > zero!() {
        aa = a;
        pp = p.sqrt();
        bb = b / pp;
    } else {
        let f = kc * kc;
        let q = (one!() - f) * (b - a * p);
        let g = one!() - p;
        let h = f - p;
        pp = (h / g).sqrt();
        aa = (a - b) / g;
        bb = -q / (g * g * pp) + aa * pp;
    }

    let mut em = one!();
    let mut f = aa;
    aa = aa + bb / pp;
    let mut g = k / pp;
    bb = two!() * (bb + f * g);
    pp = pp + g;
    g = em;
    em = em + k;
    let mut kk = k;

    while (g - k).abs() > (g * T::ca()) {
        k = two!() * kk.sqrt();
        kk = k * em;
        f = aa;
        aa = aa + bb / pp;
        g = kk / pp;
        bb = two!() * (bb + f * g);
        pp = pp + g;
        g = em;
        em = em + k;
    }

    Ok(pi_2!() * (bb + aa * em) / (em * (em + pp)))
}

// Reference: Bulirsch, “Numerical Calculation of Elliptic Integrals and Elliptic Functions.”
/// Compute [complete elliptic integral of the first kind in Bulirsch's form](https://link.springer.com/article/10.1007/bf01397975).
/// ```text
///               π/2                                                   
///              ⌠               dϑ              
/// cel1(kc)  =  ⎮  ────────────────────────────
///              ⎮      ______________________
///              ⌡   ╲╱ cos²(ϑ) + kc² sin²(ϑ)    
///             0                                                   
/// where kc ≠ 0
/// ```
///
/// Note that kc² = mc = 1 - m.
pub fn cel1<T: Float + BulirschConst>(kc: T) -> Result<T, &'static str> {
    if kc == zero!() {
        return Err("cel1: kc cannot be zero.");
    }

    let mut kc = kc.abs();
    let mut m = one!();

    loop {
        let h = m;
        m = kc + m;

        if (h - kc) > T::ca() * h {
            kc = (h * kc).sqrt();
            m = m / two!();
            continue;
        }

        break;
    }

    Ok(pi!() / m)
}

// Reference: Bulirsch, “Numerical Calculation of Elliptic Integrals and Elliptic Functions.”
/// Compute [complete elliptic integral of the second kind in Bulirsch's form](https://link.springer.com/article/10.1007/bf01397975)
/// ```text
///                     π/2                           
///                    ⌠              a + b tan²(ϑ)
/// cel2(kc, a, b)  =  ⎮  ──────────────────────────────────── ⋅ dϑ
///                    ⎮     ________________________________
///                    ⌡  ╲╱ (1 + tan²(ϑ)) (1 + kc² tan²(ϑ))    
///                    0                                                   
/// where kc ≠ 0
/// ```
///
/// Note that kc² = mc = 1 - m.
pub fn cel2<T: Float + BulirschConst>(kc: T, a: T, b: T) -> Result<T, &'static str> {
    if kc == zero!() {
        return Err("cel2: kc cannot be zero.");
    }

    let mut kc = kc.abs();
    let mut a = a;
    let mut b = b;

    let mut m = one!();
    let mut c = a;
    a = b + a;

    loop {
        b = (c * kc + b) * two!();
        c = a;
        let m0 = m;
        m = kc + m;
        a = b / m + a;

        if (m0 - kc).abs() > T::ca() * m0 {
            kc = (kc * m0).sqrt() * two!();
            continue;
        }

        break;
    }

    Ok(pi!() / four!() * a / m)
}

#[cfg(test)]
mod tests {
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
            assert_close!(ellipk, cel(kc, 1.0, 1.0, 1.0).unwrap(), 5e-12);
            assert_close!(ellipe, cel(kc, 1.0, 1.0, kc * kc).unwrap(), 1e-14);
            assert_close!(
                (ellipe - kc * kc * ellipk) / m,
                cel(kc, 1.0, 1.0, 0.0).unwrap(),
                2e-14
            );
        }

        let linsp_neg = linspace(-1.0, -1e-3, 100);
        linsp_neg.iter().for_each(|kc| test_kc(*kc));
        let linsp_pos = linspace(1e-3, 1.0, 100);
        linsp_pos.iter().for_each(|kc| test_kc(*kc));
    }

    #[test]
    fn test_cel1() {
        fn test_kc(kc: f64) {
            let k = (1.0 - kc * kc).sqrt();
            let m = k * k;
            assert_close!(ellipk(m).unwrap(), cel1(kc).unwrap(), 5e-12);
        }

        let linsp_neg = linspace(-1.0, -1e-3, 100);
        linsp_neg.iter().for_each(|kc| test_kc(*kc));
        let linsp_pos = linspace(1e-3, 1.0, 100);
        linsp_pos.iter().for_each(|kc| test_kc(*kc));
    }

    #[test]
    fn test_cel2() {
        fn test_kc(kc: f64) {
            let k = (1.0 - kc * kc).sqrt();
            let m = k * k;
            assert_close!(ellipe(m).unwrap(), cel2(kc, 1.0, kc * kc).unwrap(), 7.7e-16);
        }

        let linsp_neg = linspace(-1.0, -1e-3, 100);
        linsp_neg.iter().for_each(|kc| test_kc(*kc));
        let linsp_pos = linspace(1e-3, 1.0, 100);
        linsp_pos.iter().for_each(|kc| test_kc(*kc));
    }
}
