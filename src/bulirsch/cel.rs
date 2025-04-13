/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use num_traits::Float;

use super::BulirschConst;

// Reference: Bulirsch, “Numerical Calculation of Elliptic Integrals and Elliptic Functions III.”
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
/// Note that kc² = mc = 1 - m. The precision can be adjusted by overwriting the trait [super::BulirschConst].
///
/// # Examples
/// ```
/// use ellip::{cel, util::assert_close};
///
/// assert_close(cel(0.5, 1.0, 1.0, 1.0).unwrap(), 2.1565156474996434, 1e-15);
/// ```
pub fn cel<T: Float + BulirschConst>(kc: T, p: T, a: T, b: T) -> Result<T, &'static str> {
    if kc == zero!() {
        return Err("cel: kc cannot be zero.");
    }

    if p == zero!() {
        return Err("cel: p cannot be zero.");
    }

    let mut kc = kc.abs();
    let mut p = p;
    let mut a = a;
    let mut b = b;

    let mut e = kc;
    let mut m = one!();

    if p > zero!() {
        p = p.sqrt();
        b = b / p;
    } else {
        let mut f = kc * kc;
        let mut q = one!() - f;
        let g = one!() - p;
        f = f - p;
        q = (b - a * p) * q;
        p = (f / g).sqrt();
        a = (a - b) / g;
        b = -q / (g * g * p) + a * p;
    }

    loop {
        let f = a;
        a = b / p + a;
        let g = e / p;
        b = two!() * (f * g + b);
        p = g + p;
        let g = m;
        m = kc + m;

        if (g - kc).abs() > g * T::ca() {
            kc = two!() * e.sqrt();
            e = kc * m;
            continue;
        }

        break;
    }

    Ok(pi_2!() * (a * m + b) / (m * (m + p)))
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
/// Note that kc² = mc = 1 - m. The precision can be adjusted by overwriting the trait [super::BulirschConst].
///
/// # Examples
/// ```
/// use ellip::{cel1, util::assert_close};
///
/// assert_close(cel1(0.5).unwrap(), 2.1565156474996434, 1e-15);
/// ```
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
/// Compute [complete elliptic integral of the second kind in Bulirsch's form](https://link.springer.com/article/10.1007/bf01397975).
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
/// Note that kc² = mc = 1 - m. The precision can be adjusted by overwriting the trait [super::BulirschConst].
///
/// # Examples
/// ```
/// use ellip::{cel2, util::assert_close};
///
/// assert_close(cel2(0.5, 1.0, 1.0).unwrap(), 2.1565156474996434, 1e-15);
/// ```
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
    use itertools::iproduct;

    use super::*;
    use crate::{assert_close, ellipe, ellipk, ellippi, test_util::linspace};

    /// Test using relationship with Legendre form.
    /// Reference: https://dlmf.nist.gov/19.2#iii
    #[test]
    fn test_cel() {
        fn _test(kc: f64, p: f64) {
            let m = 1.0 - kc * kc;
            let ellipk = ellipk(m).unwrap();
            let ellipe = ellipe(m).unwrap();

            // cel precision is low for K cases
            assert_close!(ellipk, cel(kc, 1.0, 1.0, 1.0).unwrap(), 2e-12);
            assert_close!(ellipe, cel(kc, 1.0, 1.0, kc * kc).unwrap(), 1e-15);
            assert_close!(
                (ellipe - kc * kc * ellipk) / m,
                cel(kc, 1.0, 1.0, 0.0).unwrap(),
                8.5e-14
            );

            // Bulirsch, “Numerical Calculation of Elliptic Integrals and Elliptic Functions III”
            let n = 1.0 - p;
            let ellippi = ellippi(n, m).unwrap();

            assert_close!(
                (ellipk - ellipe) / m,
                cel(kc, 1.0, 0.0, 1.0).unwrap(),
                6e-12
            );
            // cel precision is very low for PI cases
            assert_close!(ellippi, cel(kc, p, 1.0, 1.0).unwrap(), 3.5e-12);
            assert_close!(
                (ellippi - ellipk) / (1.0 - p),
                cel(kc, p, 0.0, 1.0).unwrap(),
                3.5e-12
            );
        }

        let linsp_kc = [
            linspace(-1.0 + 1e-3, -1e-3, 100),
            linspace(1e-3, 1.0 - 1e-3, 100),
        ]
        .concat();
        let linsp_p = linspace(1e-3, 1.0 - 1e-3, 10);

        iproduct!(linsp_kc, linsp_p).for_each(|(kc, p)| _test(kc, p));

        // Data from Bulirsch, “Numerical Calculation of Elliptic Integrals and Elliptic Functions III”
        assert_close!(cel(1e-1, 4.1, 1.2, 1.1).unwrap(), 1.5464442694017956, 5e-16);
        assert_close!(
            cel(1e-1, -4.1, 1.2, 1.1).unwrap(),
            -6.7687378198360556e-1,
            5e-16
        );
    }

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
