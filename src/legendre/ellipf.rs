/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is translated from SciPy C++ implementation to Rust.
 */

/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     ellik.c
 *
 *     Incomplete elliptic integral of the first kind
 *
 *
 *
 * SYNOPSIS:
 *
 * double phi, m, y, ellik();
 *
 * y = ellik( phi, m );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *
 *                phi
 *                 -
 *                | |
 *                |           dt
 * F(phi | m) =   |    ------------------
 *                |                   2
 *              | |    sqrt( 1 - m sin t )
 *               -
 *                0
 *
 * of amplitude phi and modulus m, using the arithmetic -
 * geometric mean algorithm.
 *
 *
 *
 *
 * ACCURACY:
 *
 * Tested at random points with m in [0, 1] and phi as indicated.
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     -10,10       200000      7.4e-16     1.0e-16
 *
 *
 */

/*
 * Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1984, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */
/* Copyright 2014, Eric W. Moore */

use std::f64::consts::{FRAC_PI_2, PI};

use num_traits::Float;

use crate::ellipk;

/// Compute [incomplete elliptic integral of the first kind](https://dlmf.nist.gov/19.2.E4).
///
/// ```text
///              φ
///             ⌠          dθ
/// F(φ, m)  =  │  _________________
///             │     _____________
///             ⌡   \╱ 1 - m sin²θ
///            0
/// where m ≤ 1
/// ```
///
/// Note that some mathematical references use the parameter k for the function,
/// where k² = m.
pub fn ellipf<T: Float>(phi: T, m: T) -> Result<T, &'static str> {
    if m > T::one() {
        return Err("ellipf: m must be less than 1.");
    }

    if phi.is_infinite() || m.is_infinite() {
        if phi.is_infinite() {
            return Ok(T::infinity());
        }

        if m.is_infinite() {
            return Ok(T::zero());
        }

        return Err("ellipf: m or φ must be finite");
    }

    if m == T::zero() {
        return Ok(phi);
    }

    if m == T::one() {
        if phi.abs() >= T::from(FRAC_PI_2).unwrap() {
            return Ok(T::infinity());
        }
        return Ok(phi.tan().sinh().atanh());
    }

    let mut npio2 = (phi / T::from(FRAC_PI_2).unwrap()).floor();
    if (npio2.abs() % T::from(2.0).unwrap()) == T::one() {
        npio2 = npio2 + T::one();
    }

    let mut phi = phi;
    let mut k;
    if npio2 != T::zero() {
        k = ellipk(m)?;
        phi = phi - npio2 * T::from(FRAC_PI_2).unwrap();
    } else {
        k = T::zero();
    }

    let sign = if phi < T::zero() {
        phi = -phi;
        -T::one()
    } else {
        T::one()
    };

    fn done<T: Float>(val: T, sign: T, npio2: T, k: T) -> Result<T, &'static str> {
        Ok(sign * val + npio2 * k)
    }

    let a = T::one() - m;
    if a > T::one() {
        return done(ellipf_neg_m(phi, m), sign, npio2, k);
    }

    let mut b = a.sqrt();
    let mut t = phi.tan();
    if t.abs() > T::from(10.0).unwrap() {
        let mut e = T::one() / (b * t);
        if e.abs() < T::from(10.0).unwrap() {
            e = e.atan();
            if npio2 == T::zero() {
                k = ellipk(m)?;
            }
            return done(k - ellipf(e, m)?, sign, npio2, k);
        }
    }

    let mut a = T::one();
    let mut c = m.sqrt();
    let mut d = T::one();
    let mut mod_phi = T::zero();

    while (c / a).abs() > T::epsilon() {
        let temp = b / a;
        phi = phi + (t * temp).atan() + mod_phi * T::from(PI).unwrap();
        let denom = T::one() - temp * t * t;
        if denom.abs() > T::from(10.0).unwrap() * T::epsilon() {
            t = t * (T::one() + temp) / denom;
            mod_phi = T::from(
                ((phi + T::from(FRAC_PI_2).unwrap()) / T::from(PI).unwrap())
                    .to_i32()
                    .unwrap(),
            )
            .unwrap();
        } else {
            t = phi.tan();
            mod_phi = T::from(
                ((phi - t.atan()) / T::from(PI).unwrap())
                    .floor()
                    .to_i32()
                    .unwrap(),
            )
            .unwrap();
        }
        c = (a - b) / T::from(2.0).unwrap();
        let temp = (a * b).sqrt();
        a = (a + b) / T::from(2.0).unwrap();
        b = temp;
        d = d * T::from(2.0).unwrap();
    }

    done(
        (t.atan() + mod_phi * T::from(PI).unwrap()) / (d * a),
        sign,
        npio2,
        k,
    )
}

/// Compute elliptic integral of the first kind for m<0.
#[inline]
fn ellipf_neg_m<T: Float>(phi: T, m: T) -> T {
    let mpp = m * phi * phi;

    if -mpp < T::from(1e-6).unwrap() && phi < -m {
        return phi
            + (-mpp * phi * phi / T::from(30.0).unwrap()
                + T::from(3.0).unwrap() * mpp * mpp / T::from(40.0).unwrap()
                + mpp / T::from(6.0).unwrap())
                * phi;
    }

    if -mpp > T::from(4e7).unwrap() {
        let sm = (-m).sqrt();
        let sp = phi.sin();
        let cp = phi.cos();

        let a = (T::from(4.0).unwrap() * sp * sm / (T::one() + cp)).ln();
        let b = -(T::one() + cp / (sp * sp) - a) / (T::from(4.0).unwrap() * m);
        return (a + b) / sm;
    }

    let (scale, x, y, z) = if phi > T::from(1e-153).unwrap() && m > T::from(-1e305).unwrap() {
        let s = phi.sin();
        let phi_tan = phi.tan();
        let csc2 = T::one() / (s * s);
        (T::one(), T::one() / (phi_tan * phi_tan), csc2 - m, csc2)
    } else {
        (phi, T::one(), T::one() - m * phi * phi, T::one())
    };

    if x == y && x == z {
        return scale / x.sqrt();
    }

    let a0 = (x + y + z) / T::from(3.0).unwrap();
    let mut a = a0;
    let mut x1 = x;
    let mut y1 = y;
    let mut z1 = z;
    let mut q = T::from(400.0).unwrap() * (a0 - x).abs().max((a0 - y).abs()).max((a0 - z).abs());
    let mut n = 0;

    while q > a.abs() && n <= 100 {
        let sx = x1.sqrt();
        let sy = y1.sqrt();
        let sz = z1.sqrt();
        let lam = sx * sy + sx * sz + sy * sz;
        x1 = (x1 + lam) / T::from(4.0).unwrap();
        y1 = (y1 + lam) / T::from(4.0).unwrap();
        z1 = (z1 + lam) / T::from(4.0).unwrap();
        a = (x1 + y1 + z1) / T::from(3.0).unwrap();
        n += 1;
        q = q / T::from(4.0).unwrap();
    }

    let two_to_2n = T::from(1u32 << (2 * n)).unwrap();
    let x = (a0 - x) / a / two_to_2n;
    let y = (a0 - y) / a / two_to_2n;
    let z = -(x + y);

    let e2 = x * y - z * z;
    let e3 = x * y * z;

    scale
        * (T::one() - e2 / T::from(10.0).unwrap()
            + e3 / T::from(14.0).unwrap()
            + e2 * e2 / T::from(24.0).unwrap()
            - T::from(3.0).unwrap() * e2 * e3 / T::from(44.0).unwrap())
        / a.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compare_test_data;

    fn ellipf_k(inp: &[f64]) -> f64 {
        ellipf(inp[0], inp[1] * inp[1]).unwrap()
    }

    #[test]
    fn test_ellipf() {
        compare_test_data!("./tests/data/boost/ellint_f_data.txt", ellipf_k, 4.9e-16);
    }
}
