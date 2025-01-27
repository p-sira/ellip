/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is translated from SciPy C++ implementation to Rust.
 */

/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     ellie.c
 *
 *     Incomplete elliptic integral of the second kind
 *
 *
 *
 * SYNOPSIS:
 *
 * double phi, m, y, ellie();
 *
 * y = ellie( phi, m );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *                 phi
 *                  -
 *                 | |
 *                 |                   2
 * E(phi_\m)  =    |    sqrt( 1 - m sin t ) dt
 *                 |
 *               | |
 *                -
 *                 0
 *
 * of amplitude phi and modulus m, using the arithmetic -
 * geometric mean algorithm.
 *
 *
 *
 * ACCURACY:
 *
 * Tested at random arguments with phi in [-10, 10] and m in
 * [0, 1].
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     -10,10      150000       3.3e-15     1.4e-16
 */

/*
 * Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1984, 1987, 1993 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */
/* Copyright 2014, Eric W. Moore */
use std::f64::consts::{FRAC_PI_2, PI};

use crate::{ellipe, ellipk};

/// Compute incomplete elliptic integral of the second kind.
/// ```text
///              φ
///             ⌠     _______________
/// E(φ, m)  =  │   \╱ 1 - m sin²(t)  dt
///             ⌡
///            0
/// where m ≤ 1
/// ```
pub fn ellipeinc(phi: f64, m: f64) -> Result<f64, &'static str> {
    if m > 1.0 {
        return Err("ellipeinc: m must satisfy: m ≤ 1.");
    }

    if phi.is_infinite() {
        return Ok(phi);
    }

    if m.is_infinite() {
        return Ok(f64::NEG_INFINITY);
    }

    if m == 0.0 {
        return Ok(phi);
    }

    let mut lphi = phi;
    let mut npio2 = (lphi / FRAC_PI_2).floor();

    if (npio2.abs() % 2.0) == 1.0 {
        npio2 += 1.0;
    }

    lphi -= npio2 * FRAC_PI_2;

    let sign = if lphi < 0.0 {
        lphi = -lphi;
        -1.0
    } else {
        1.0
    };

    let a = 1.0 - m;
    let e = ellipe(m)?;

    fn done(val: f64, sign: f64, npio2: f64, e: f64) -> Result<f64, &'static str> {
        Ok(sign * val + npio2 * e)
    }

    if a == 0.0 {
        return done(lphi.sin(), sign, npio2, e);
    }

    if a > 1.0 {
        return done(ellipeinc_neg_m(lphi, m), sign, npio2, e);
    }

    if lphi < 0.135 {
        let m11 = (((((-7.0 / 2816.0) * m + (5.0 / 1056.0)) * m - (7.0 / 2640.0)) * m
            + (17.0 / 41580.0))
            * m
            - (1.0 / 155925.0))
            * m;
        let m9 =
            ((((-5.0 / 1152.0) * m + (1.0 / 144.0)) * m - (1.0 / 360.0)) * m + (1.0 / 5670.0)) * m;
        let m7 = ((-m / 112.0 + (1.0 / 84.0)) * m - (1.0 / 315.0)) * m;
        let m5 = (-m / 40.0 + (1.0 / 30.0)) * m;
        let m3 = -m / 6.0;
        let p2 = lphi * lphi;

        return done(
            ((((m11 * p2 + m9) * p2 + m7) * p2 + m5) * p2 + m3) * p2 * lphi + lphi,
            sign,
            npio2,
            e,
        );
    }

    let mut t = lphi.tan();
    let mut b = a.sqrt();

    if t.abs() > 10.0 {
        let e = 1.0 / (b * t);
        if e.abs() < 10.0 {
            let e = e.atan();
            return done(
                e + m * lphi.sin() * e.sin() - ellipeinc(e, m)?,
                sign,
                npio2,
                e,
            );
        }
    }

    let mut c = m.sqrt();
    let mut a = 1.0;
    let mut d = 1.0;
    let mut ee = 0.0;
    let mut mod_phi = 0.0;

    while (c / a).abs() > f64::EPSILON {
        let temp = b / a;
        lphi += (t * temp).atan() + mod_phi * PI;
        let denom = 1.0 - temp * t * t;

        if denom.abs() > 10.0 * f64::EPSILON {
            t *= (1.0 + temp) / denom;
            mod_phi = ((lphi + FRAC_PI_2) / PI) as i32 as f64;
        } else {
            t = lphi.tan();
            mod_phi = ((lphi - t.atan()) / PI).floor() as i32 as f64;
        }

        c = (a - b) / 2.0;
        let temp = (a * b).sqrt();
        a = (a + b) / 2.0;
        b = temp;
        d += d;
        ee += c * lphi.sin();
    }

    done(
        e / ellipk(m)? * (t.atan() + mod_phi * PI) / (d * a) + ee,
        sign,
        npio2,
        e,
    )
}

/// Compute elliptic integral of the second kind for m<0.
fn ellipeinc_neg_m(phi: f64, m: f64) -> f64 {
    let mpp = m * phi * phi;

    if -mpp < 1e-6 && phi < -m {
        return phi + (mpp * phi * phi / 30.0 - mpp * mpp / 40.0 - mpp / 6.0) * phi;
    }

    if -mpp > 1e6 {
        let sm = (-m).sqrt();
        let sp = phi.sin();
        let cp = phi.cos();

        let a = 1.0 - cp;
        let b1 = (4.0 * sp * sm / (1.0 + cp)).ln();
        let b = -(0.5 + b1) / 2.0 / m;
        let c = (0.75 + cp / sp / sp - b1) / 16.0 / m / m;
        return (a + b + c) * sm;
    }

    let scalef: f64;
    let scaled: f64;
    let x: f64;
    let y: f64;
    let z: f64;
    if phi > 1e-153 && m > -1e200 {
        let s = phi.sin();
        let csc2 = 1.0 / (s * s);
        scalef = 1.0;
        scaled = m / 3.0;
        let phi_tan = phi.tan();
        x = 1.0 / (phi_tan * phi_tan);
        y = csc2 - m;
        z = csc2;
    } else {
        scalef = phi;
        scaled = mpp * phi / 3.0;
        x = 1.0;
        y = 1.0 - mpp;
        z = 1.0;
    }

    if x == y && x == z {
        return (scalef + scaled / x) / x.sqrt();
    }

    let a0f: f64 = (x + y + z) / 3.0;
    let a0d: f64 = (x + y + 3.0 * z) / 5.0;
    let mut af = a0f;
    let mut ad = a0d;

    let mut seriesd = 0.0;
    let mut seriesn = 1.0;

    let mut q = 400.0 * (a0f - x).abs().max((a0f - y).abs().max((a0f - z).abs()));

    let mut x1 = x;
    let mut y1 = y;
    let mut z1 = z;
    let mut n: i32 = 0;
    while q > af.abs() && q > ad.abs() && n <= 100 {
        let sx = x1.sqrt();
        let sy = y1.sqrt();
        let sz = z1.sqrt();
        let lam = sx * sy + sx * sz + sy * sz;
        seriesd += seriesn / (sz * (z1 + lam));
        x1 = (x1 + lam) / 4.0;
        y1 = (y1 + lam) / 4.0;
        z1 = (z1 + lam) / 4.0;
        af = (x1 + y1 + z1) / 3.0;
        ad = (ad + lam) / 4.0;
        n += 1;
        q /= 4.0;
        seriesn /= 4.0;
    }

    let two_to_2n = (1 << (2 * n)) as f64;
    let xf: f64 = (a0f - x) / af / two_to_2n;
    let yf: f64 = (a0f - y) / af / two_to_2n;
    let zf: f64 = -(xf + yf);

    let e2f: f64 = xf * yf - zf * zf;
    let e3f: f64 = xf * yf * zf;

    let mut ret = scalef
        * (1.0 - e2f / 10.0 + e3f / 14.0 + e2f * e2f / 24.0 - 3.0 * e2f * e3f / 44.0)
        / af.sqrt();

    let xd: f64 = (a0d - x) / ad / two_to_2n;
    let yd: f64 = (a0d - y) / ad / two_to_2n;
    let zd: f64 = -(xd + yd) / 3.0;

    let e2d: f64 = xd * yd - 6.0 * zd * zd;
    let e3d: f64 = (3.0 * xd * yd - 8.0 * zd * zd) * zd;
    let e4d: f64 = 3.0 * (xd * yd - zd * zd) * zd * zd;
    let e5d: f64 = xd * yd * zd * zd * zd;

    ret -= scaled
        * (1.0 - 3.0 * e2d / 14.0 + e3d / 6.0 + 9.0 * e2d * e2d / 88.0
            - 3.0 * e4d / 22.0
            - 9.0 * e2d * e3d / 52.0
            + 3.0 * e5d / 26.0)
        / two_to_2n
        / ad
        / ad.sqrt();

    ret -= 3.0 * scaled * seriesd;
    ret
}

#[cfg(test)]
mod test {
    use crate::compare_test_data;

    use super::*;

    fn ellipeinc_k(inp: &[f64]) -> f64 {
        ellipeinc(inp[0], inp[1] * inp[1]).unwrap()
    }

    #[test]
    fn test_ellipeinc() {
        compare_test_data!("./tests/data/boost/ellint_e2_data.txt", ellipeinc_k, 1e-14);
    }
}
