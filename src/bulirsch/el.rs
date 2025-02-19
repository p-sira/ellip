/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use num_traits::Float;

use super::BulirschConst;

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
    return Ok(if x < zero!() { -e } else { e });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{assert_close, ellipf, test_util::linspace};

    #[test]
    fn test_el1() {
        fn _test(x: f64, kc: f64) {
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
            .for_each(|(&x, &kc)| _test(x, kc));
        let linsp_pos = linspace(1e-3, 1.0, 50);
        x.iter()
            .zip(linsp_pos.iter())
            .for_each(|(&x, &kc)| _test(x, kc));
    }
}
