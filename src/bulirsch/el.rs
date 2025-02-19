/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use num_traits::Float;

use super::{cel1, cel2, BulirschConst};

// Reference: Bulirsch, “Numerical Calculation of Elliptic Integrals and Elliptic Functions.”
/// Compute [incomplete elliptic integral of the first kind in Bulirsch's form](https://dlmf.nist.gov/19.2.E11_5)
/// ```text
///                 arctan(x)                                                   
///                ⌠                             
///                |               dϑ
/// el1(x, kc)  =  ⎮  ────────────────────────────
///                ⎮      ______________________
///                ⌡   ╲╱ cos²(ϑ) + kc² sin²(ϑ)    
///               0                                                   
/// where kc ≠ 0
/// ```
///
/// Note that x = tan φ and kc² = mc = 1 - m.
pub fn el1<T: Float + BulirschConst>(x: T, kc: T) -> Result<T, &'static str> {
    if x == zero!() {
        return Ok(zero!());
    }

    if kc == zero!() {
        return Err("el1: kc cannot be zero.");
    }

    // phi = π/2
    if x == inf!() {
        return cel1(kc);
    }

    let mut y = (one!() / x).abs();
    let mut kc = kc.abs();
    let mut m = one!();
    let mut l = 0;

    loop {
        let e = m * kc;
        let g = m;
        m = kc + m;
        y = -e / y + y;

        if y == zero!() {
            y = e.sqrt() * T::cb();
        }

        if (g - kc).abs() > T::ca() * g {
            kc = e.sqrt() * two!();
            l *= 2;
            if y < zero!() {
                l += 1;
            }
            continue;
        }

        break;
    }

    if y < zero!() {
        l += 1;
    }

    let e = ((m / y).atan() + pi!() * num!(l)) / m;
    Ok(if x < zero!() { -e } else { e })
}

// Reference: Bulirsch, “Numerical Calculation of Elliptic Integrals and Elliptic Functions.”
/// Compute [incomplete elliptic integral of the second kind in Bulirsch's form](https://dlmf.nist.gov/19.2.E12)
/// ```text
///                       arctan(x)                                                   
///                      ⌠                             
///                      |              a + b tan²(ϑ)
/// el1(x, kc, a, b)  =  ⎮  ──────────────────────────────────── ⋅ dϑ
///                      ⎮     ________________________________
///                      ⌡  ╲╱ (1 + tan²(ϑ)) (1 + kc² tan²(ϑ))    
///                     0                                                   
/// where kc ≠ 0
/// ```
///
/// Note that x = tan φ and kc² = mc = 1 - m.
pub fn el2<T: Float + BulirschConst>(x: T, kc: T, a: T, b: T) -> Result<T, &'static str> {
    if x == zero!() {
        return Ok(zero!());
    }

    if kc == zero!() {
        return Err("el2: kc cannot be zero.");
    }

    // phi = π/2
    if x == inf!() {
        return cel2(kc, a, b);
    }

    let mut b = b;
    let mut c = x * x;
    let mut d = one!() + c;
    let mut p = ((one!() + kc * kc * c) / d).sqrt();

    d = x / d;
    c = d / (p * two!());
    let z = a - b;
    let mut i = a;
    let mut a = (b + a) / two!();
    let mut y = (one!() / x).abs();
    let mut f = zero!();
    let mut l = 0;
    let mut m = one!();
    let mut kc = kc.abs();

    loop {
        b = i * kc + b;
        let e = m * kc;
        let mut g = e / p;
        d = f * g + d;
        f = c;
        i = a;
        p = g + p;
        c = (d / p + c) / two!();
        g = m;
        m = kc + m;
        a = (b / m + a) / two!();
        y = -e / y + y;

        if y == zero!() {
            y = e.sqrt() * T::cb();
        }

        if (g - kc).abs() > T::ca() * g {
            kc = e.sqrt() * two!();
            l *= 2;
            if y < zero!() {
                l += 1;
            }
            continue;
        }

        break;
    }

    if y < zero!() {
        l += 1;
    }

    let mut e = ((m / y).atan() + pi!() * num!(l)) * a / m;

    if x < zero!() {
        e = -e;
    }

    Ok(e + c * z)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{assert_close, ellipdinc, ellipeinc, ellipf, test_util::linspace};

    #[test]
    fn test_el1() {
        fn test_special_cases(x: f64, kc: f64) {
            let k = (1.0 - kc * kc).sqrt();
            let m = k * k;
            let phi = x.atan();
            let f = ellipf(phi, m).unwrap();

            assert_close!(f, el1(x, kc).unwrap(), 5e-14);
        }

        let x = linspace(0.0, 100.0, 50);
        let linsp_neg = linspace(-1.0, -1e-3, 50);
        x.iter()
            .zip(linsp_neg.iter())
            .for_each(|(&x, &kc)| test_special_cases(x, kc));
        let linsp_pos = linspace(1e-3, 1.0, 50);
        x.iter()
            .zip(linsp_pos.iter())
            .for_each(|(&x, &kc)| test_special_cases(x, kc));

        // Test computed values from the reference
        // Bulirsch, “Numerical Calculation of Elliptic Integrals and Elliptic Functions.”
        fn test_reference(x: f64, expected: f64) {
            // The reference corrects to 10 decimals.
            let significants = 1e10;

            assert_eq!(
                expected,
                (el1(x, 1e-11).unwrap() * significants).round() / significants
            );
        }

        test_reference(1e5, 12.2060726456);
        test_reference(1e6, 14.5086577385);
        test_reference(1e7, 16.8112428290);
        test_reference(1e8, 19.1138276745);
        test_reference(1e9, 21.4163880184);
        test_reference(1e10, 23.7165074338);
        test_reference(1e11, 25.8333567970);
        test_reference(1e12, 26.6148963052);
        test_reference(1e23, 26.7147303841);
        test_reference(f64::infinity(), 26.7147303841);
    }

    #[test]
    fn test_el2() {
        fn test_special_cases(x: f64, kc: f64) {
            let k = (1.0 - kc * kc).sqrt();
            let m = k * k;
            let phi = x.atan();
            let f = ellipf(phi, m).unwrap();
            let e = ellipeinc(phi, m).unwrap();
            let d = ellipdinc(phi, m).unwrap();

            assert_close!(f, el2(x, kc, 1.0, 1.0).unwrap(), 5e-14);
            assert_close!(e, el2(x, kc, 1.0, kc * kc).unwrap(), 7e-16);
            assert_close!(d, el2(x, kc, 0.0, 1.0).unwrap(), 6e-14);
        }

        let x = linspace(0.0, 100.0, 50);
        let linsp_neg = linspace(-1.0, -1e-3, 50);
        x.iter()
            .zip(linsp_neg.iter())
            .for_each(|(&x, &kc)| test_special_cases(x, kc));
        let linsp_pos = linspace(1e-3, 1.0, 50);
        x.iter()
            .zip(linsp_pos.iter())
            .for_each(|(&x, &kc)| test_special_cases(x, kc));
    }
}
