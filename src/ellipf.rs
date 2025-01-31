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

use crate::ellipk;

/// Compute [incomplete elliptic integral of the first kind](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipkinc.html).
///
/// ```text
///              φ
///             ⌠          dt
/// F(φ, m)  =  │  ___________________
///             │     _______________
///             ⌡   \╱ 1 - m sin²(t)
///            0
/// where m ≤ 1
/// ```
///
/// Note that some mathematical references use the parameter k for the function,
/// where k² = m.
pub fn ellipf(phi: f64, m: f64) -> Result<f64, &'static str> {
    if m > 1.0 {
        return Err("ellipf: m must satisfy: m ≤ 1.");
    }

    if phi.is_infinite() && m.is_infinite() {
        return Err("ellipf: m or φ must be finite");
    }

    if phi.is_infinite() {
        return Ok(f64::INFINITY);
    }

    if m.is_infinite() {
        return Ok(0.0);
    }

    if m == 0.0 {
        return Ok(phi);
    }

    if m == 1.0 {
        if phi.abs() >= FRAC_PI_2 {
            return Ok(f64::INFINITY);
        }
        return Ok(phi.tan().asinh());
    }

    let mut npio2 = (phi / FRAC_PI_2).floor();
    if (npio2.abs() % 2.0) == 1.0 {
        npio2 += 1.0;
    }

    let mut phi = phi;
    let mut k;
    if npio2 != 0.0 {
        k = ellipk(m)?;
        phi -= npio2 * FRAC_PI_2;
    } else {
        k = 0.0;
    }

    let sign = if phi < 0.0 {
        phi = -phi;
        -1.0
    } else {
        1.0
    };

    fn done(val: f64, sign: f64, npio2: f64, k: f64) -> Result<f64, &'static str> {
        Ok(sign * val + npio2 * k)
    }

    let a = 1.0 - m;
    if a > 1.0 {
        return done(ellipf_neg_m(phi, m), sign, npio2, k);
    }

    let mut b = a.sqrt();
    let mut t = phi.tan();
    if t.abs() > 10.0 {
        let mut e = 1.0 / (b * t);
        if e.abs() < 10.0 {
            e = e.atan();
            if npio2 == 0.0 {
                k = ellipk(m)?;
            }
            return done(k - ellipf(e, m)?, sign, npio2, k);
        }
    }

    let mut a = 1.0;
    let mut c = m.sqrt();
    let mut d = 1.0;
    let mut mod_phi = 0.0;

    while (c / a).abs() > f64::EPSILON {
        let temp = b / a;
        phi += (t * temp).atan() + mod_phi * PI;
        let denom = 1.0 - temp * t * t;
        if denom.abs() > 10.0 * f64::EPSILON {
            t *= (1.0 + temp) / denom;
            mod_phi = ((phi + FRAC_PI_2) / PI) as i32 as f64;
        } else {
            t = phi.tan();
            mod_phi = (((phi - t.atan()) / PI).floor()) as i32 as f64;
        }
        c = (a - b) / 2.0;
        let temp = (a * b).sqrt();
        a = (a + b) / 2.0;
        b = temp;
        d *= 2.0;
    }

    done((t.atan() + mod_phi * PI) / (d * a), sign, npio2, k)
}

/// Compute elliptic integral of the first kind for m<0.
#[inline]
fn ellipf_neg_m(phi: f64, m: f64) -> f64 {
    let mpp = m * phi * phi;

    if -mpp < 1e-6 && phi < -m {
        return phi + (-mpp * phi * phi / 30.0 + 3.0 * mpp * mpp / 40.0 + mpp / 6.0) * phi;
    }

    if -mpp > 4e7 {
        let sm = (-m).sqrt();
        let sp = phi.sin();
        let cp = phi.cos();

        let a = (4.0 * sp * sm / (1.0 + cp)).ln();
        let b = -(1.0 + cp / (sp * sp) - a) / (4.0 * m);
        return (a + b) / sm;
    }

    let (scale, x, y, z) = if phi > 1e-153 && m > -1e305 {
        let s = phi.sin();
        let phi_tan = phi.tan();
        let csc2 = 1.0 / (s * s);
        (1.0, 1.0 / (phi_tan * phi_tan), csc2 - m, csc2)
    } else {
        (phi, 1.0, 1.0 - m * phi * phi, 1.0)
    };

    if x == y && x == z {
        return scale / x.sqrt();
    }

    let a0 = (x + y + z) / 3.0;
    let mut a = a0;
    let mut x1 = x;
    let mut y1 = y;
    let mut z1 = z;
    let mut q = 400.0 * f64::max((a0 - x).abs(), f64::max((a0 - y).abs(), (a0 - z).abs()));
    let mut n = 0;

    while q > a.abs() && n <= 100 {
        let sx = x1.sqrt();
        let sy = y1.sqrt();
        let sz = z1.sqrt();
        let lam = sx * sy + sx * sz + sy * sz;
        x1 = (x1 + lam) / 4.0;
        y1 = (y1 + lam) / 4.0;
        z1 = (z1 + lam) / 4.0;
        a = (x1 + y1 + z1) / 3.0;
        n += 1;
        q /= 4.0;
    }

    let two_to_2n = (1 << (2 * n)) as f64;
    let x = (a0 - x) / a / two_to_2n;
    let y = (a0 - y) / a / two_to_2n;
    let z = -(x + y);

    let e2 = x * y - z * z;
    let e3 = x * y * z;

    scale * (1.0 - e2 / 10.0 + e3 / 14.0 + e2 * e2 / 24.0 - 3.0 * e2 * e3 / 44.0) / a.sqrt()
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::compare_test_data;

    fn ellipf_k(inp: &[f64]) -> f64 {
        ellipf(inp[0], inp[1] * inp[1]).unwrap()
    }

    #[test]
    fn test_ellipf() {
        compare_test_data!("./tests/data/boost/ellint_f_data.txt", ellipf_k, 1e-14);
    }
}
