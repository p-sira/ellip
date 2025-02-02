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

use num_traits::Float;

use crate::{ellipe, ellipk};

/// Compute [incomplete elliptic integral of the second kind](https://dlmf.nist.gov/19.2.E5).
/// ```text
///              φ
///             ⌠   _____________
/// E(φ, m)  =  │ \╱ 1 - m sin²θ  dθ
///             ⌡
///            0
/// where m ≤ 1
/// ```
///
/// Note that some mathematical references use the parameter k for the function,
/// where k² = m.
pub fn ellipeinc<T: Float>(phi: T, m: T) -> Result<T, &'static str> {
    if m > T::one() {
        return Err("ellipeinc: m must be less than 1.");
    }

    if phi.is_infinite() {
        return Ok(phi);
    }

    if m.is_infinite() {
        return Ok(T::neg_infinity());
    }

    if m == T::zero() {
        return Ok(phi);
    }

    let mut lphi = phi;
    let mut npio2 = (lphi / T::from(FRAC_PI_2).unwrap()).floor();

    if (npio2.abs() % T::from(2.0).unwrap()) == T::one() {
        npio2 = npio2 + T::one();
    }

    lphi = lphi - npio2 * T::from(FRAC_PI_2).unwrap();

    let sign = if lphi < T::zero() {
        lphi = -lphi;
        -T::one()
    } else {
        T::one()
    };

    let a = T::one() - m;
    let e = ellipe(m)?;

    fn done<T: Float>(val: T, sign: T, npio2: T, e: T) -> Result<T, &'static str> {
        Ok(sign * val + npio2 * e)
    }

    if a == T::zero() {
        return done(lphi.sin(), sign, npio2, e);
    }

    if a > T::one() {
        return done(ellipeinc_neg_m(lphi, m), sign, npio2, e);
    }

    if lphi < T::from(0.135).unwrap() {
        let m11 = (((((-T::from(7.0).unwrap() / T::from(2816.0).unwrap()) * m
            + (T::from(5.0).unwrap() / T::from(1056.0).unwrap()))
            * m
            - (T::from(7.0).unwrap() / T::from(2640.0).unwrap()))
            * m
            + (T::from(17.0).unwrap() / T::from(41580.0).unwrap()))
            * m
            - (T::from(1.0).unwrap() / T::from(155925.0).unwrap()))
            * m;
        let m9 = ((((-T::from(5.0).unwrap() / T::from(1152.0).unwrap()) * m
            + (T::from(1.0).unwrap() / T::from(144.0).unwrap()))
            * m
            - (T::from(1.0).unwrap() / T::from(360.0).unwrap()))
            * m
            + (T::from(1.0).unwrap() / T::from(5670.0).unwrap()))
            * m;
        let m7 =
            ((-m / T::from(112.0).unwrap() + (T::from(1.0).unwrap() / T::from(84.0).unwrap())) * m
                - (T::from(1.0).unwrap() / T::from(315.0).unwrap()))
                * m;
        let m5 =
            (-m / T::from(40.0).unwrap() + (T::from(1.0).unwrap() / T::from(30.0).unwrap())) * m;
        let m3 = -m / T::from(6.0).unwrap();
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

    if t.abs() > T::from(10.0).unwrap() {
        let e = T::one() / (b * t);
        if e.abs() < T::from(10.0).unwrap() {
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
    let mut a = T::one();
    let mut d = T::one();
    let mut ee = T::zero();
    let mut mod_phi = 0;

    while (c / a).abs() > T::epsilon() {
        let temp = b / a;
        lphi = lphi + ((t * temp).atan() + T::from(mod_phi).unwrap() * T::from(PI).unwrap());
        let denom = T::one() - temp * t * t;

        if denom.abs() > T::from(10.0).unwrap() * T::epsilon() {
            t = t * (T::one() + temp) / denom;
            mod_phi = ((lphi + T::from(FRAC_PI_2).unwrap()) / T::from(PI).unwrap())
                .to_i32()
                .unwrap();
        } else {
            t = lphi.tan();
            mod_phi = ((lphi - t.atan()) / T::from(PI).unwrap())
                .floor()
                .to_i32()
                .unwrap();
        }

        c = (a - b) / T::from(2.0).unwrap();
        let temp = (a * b).sqrt();
        a = (a + b) / T::from(2.0).unwrap();
        b = temp;
        d = d + d;
        ee = ee + c * lphi.sin();
    }

    done(
        e / ellipk(m)? * (t.atan() + T::from(mod_phi).unwrap() * T::from(PI).unwrap()) / (d * a)
            + ee,
        sign,
        npio2,
        e,
    )
}

/// Compute elliptic integral of the second kind for m<0.
#[inline]
fn ellipeinc_neg_m<T: Float>(phi: T, m: T) -> T {
    let mpp = m * phi * phi;

    if -mpp < T::from(1e-6).unwrap() && phi < -m {
        return phi
            + (mpp * phi * phi / T::from(30.0).unwrap()
                - mpp * mpp / T::from(40.0).unwrap()
                - mpp / T::from(6.0).unwrap())
                * phi;
    }

    if -mpp > T::from(1e6).unwrap() {
        let sm = (-m).sqrt();
        let sp = phi.sin();
        let cp = phi.cos();

        let a = T::one() - cp;
        let b1 = (T::from(4.0).unwrap() * sp * sm / (T::one() + cp)).ln();
        let b = -(T::from(0.5).unwrap() + b1) / T::from(2.0).unwrap() / m;
        let c = (T::from(0.75).unwrap() + cp / sp / sp - b1) / T::from(16.0).unwrap() / m / m;
        return (a + b + c) * sm;
    }

    let scalef: T;
    let scaled: T;
    let x: T;
    let y: T;
    let z: T;
    if phi > T::from(1e-153).unwrap() && m > T::from(-1e200).unwrap() {
        let s = phi.sin();
        let csc2 = T::one() / (s * s);
        scalef = T::one();
        scaled = m / T::from(3.0).unwrap();
        let phi_tan = phi.tan();
        x = T::one() / (phi_tan * phi_tan);
        y = csc2 - m;
        z = csc2;
    } else {
        scalef = phi;
        scaled = mpp * phi / T::from(3.0).unwrap();
        x = T::one();
        y = T::one() - mpp;
        z = T::one();
    }

    if x == y && x == z {
        return (scalef + scaled / x) / x.sqrt();
    }

    let a0f = (x + y + z) / T::from(3.0).unwrap();
    let a0d = (x + y + T::from(3.0).unwrap() * z) / T::from(5.0).unwrap();
    let mut af = a0f;
    let mut ad = a0d;

    let mut seriesd = T::zero();
    let mut seriesn = T::one();

    let mut q = T::from(400.0).unwrap() * (a0f - x).abs().max((a0f - y).abs().max((a0f - z).abs()));

    let mut x1 = x;
    let mut y1 = y;
    let mut z1 = z;
    let mut n: i32 = 0;
    while q > af.abs() && q > ad.abs() && n <= 100 {
        let sx = x1.sqrt();
        let sy = y1.sqrt();
        let sz = z1.sqrt();
        let lam = sx * sy + sx * sz + sy * sz;
        seriesd = seriesd + seriesn / (sz * (z1 + lam));
        x1 = (x1 + lam) / T::from(4.0).unwrap();
        y1 = (y1 + lam) / T::from(4.0).unwrap();
        z1 = (z1 + lam) / T::from(4.0).unwrap();
        af = (x1 + y1 + z1) / T::from(3.0).unwrap();
        ad = (ad + lam) / T::from(4.0).unwrap();
        n += 1;
        q = q / T::from(4.0).unwrap();
        seriesn = seriesn / T::from(4.0).unwrap();
    }

    let two_to_2n = T::from(1 << (2 * n)).unwrap();
    let xf = (a0f - x) / af / two_to_2n;
    let yf = (a0f - y) / af / two_to_2n;
    let zf = -(xf + yf);

    let e2f = xf * yf - zf * zf;
    let e3f = xf * yf * zf;

    let mut ret = scalef
        * (T::one() - e2f / T::from(10.0).unwrap()
            + e3f / T::from(14.0).unwrap()
            + e2f * e2f / T::from(24.0).unwrap()
            - T::from(3.0).unwrap() * e2f * e3f / T::from(44.0).unwrap())
        / af.sqrt();

    let xd = (a0d - x) / ad / two_to_2n;
    let yd = (a0d - y) / ad / two_to_2n;
    let zd = -(xd + yd) / T::from(3.0).unwrap();

    let e2d = xd * yd - T::from(6.0).unwrap() * zd * zd;
    let e3d = (T::from(3.0).unwrap() * xd * yd - T::from(8.0).unwrap() * zd * zd) * zd;
    let e4d = T::from(3.0).unwrap() * (xd * yd - zd * zd) * zd * zd;
    let e5d = xd * yd * zd * zd * zd;

    ret = ret
        - scaled
            * (T::one() - T::from(3.0).unwrap() * e2d / T::from(14.0).unwrap()
                + e3d / T::from(6.0).unwrap()
                + T::from(9.0).unwrap() * e2d * e2d / T::from(88.0).unwrap()
                - T::from(3.0).unwrap() * e4d / T::from(22.0).unwrap()
                - T::from(9.0).unwrap() * e2d * e3d / T::from(52.0).unwrap()
                + T::from(3.0).unwrap() * e5d / T::from(26.0).unwrap())
            / two_to_2n
            / ad
            / ad.sqrt();

    ret - T::from(3.0).unwrap() * scaled * seriesd
}

#[cfg(test)]
mod test {
    use crate::compare_test_data;

    use super::*;

    fn ellipeinc_k<T: Float>(inp: &[T]) -> T {
        ellipeinc(inp[0], inp[1] * inp[1]).unwrap()
    }

    #[test]
    fn test_ellipeinc() {
        compare_test_data!(
            "./tests/data/boost/ellint_e2_data.txt",
            ellipeinc_k::<f64>,
            1.1e-15
        );
    }
}
