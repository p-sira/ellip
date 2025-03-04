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
    if m > one!() {
        return Err("ellipeinc: m must be less than 1.");
    }

    if phi.is_infinite() {
        return Ok(phi);
    }

    if m.is_infinite() {
        return Ok(neg_inf!());
    }

    if m == zero!() {
        return Ok(phi);
    }

    let mut lphi = phi;
    let mut npio2 = (lphi / pi_2!()).floor();

    if (npio2.abs() % two!()) == one!() {
        npio2 = npio2 + one!();
    }

    lphi = lphi - npio2 * pi_2!();

    let sign = if lphi < zero!() {
        lphi = -lphi;
        -one!()
    } else {
        one!()
    };

    let a = one!() - m;
    let ellipe = ellipe(m)?;

    fn done<T: Float>(val: T, sign: T, npio2: T, ellipe: T) -> Result<T, &'static str> {
        Ok(sign * val + npio2 * ellipe)
    }

    if a == zero!() {
        return done(lphi.sin(), sign, npio2, ellipe);
    }

    if a > one!() {
        return done(ellipeinc_neg_m(lphi, m), sign, npio2, ellipe);
    }

    if lphi < num!(0.135) {
        let m11 = (((((-num!(7.0) / num!(2816.0)) * m + (five!() / num!(1056.0))) * m
            - (num!(7.0) / num!(2640.0)))
            * m
            + (num!(17.0) / num!(41580.0)))
            * m
            - (num!(1.0) / num!(155925.0)))
            * m;
        let m9 = ((((-five!() / num!(1152.0)) * m + (num!(1.0) / num!(144.0))) * m
            - (num!(1.0) / num!(360.0)))
            * m
            + (num!(1.0) / num!(5670.0)))
            * m;
        let m7 =
            ((-m / num!(112.0) + (num!(1.0) / num!(84.0))) * m - (num!(1.0) / num!(315.0))) * m;
        let m5 = (-m / num!(40.0) + (num!(1.0) / num!(30.0))) * m;
        let m3 = -m / six!();
        let p2 = lphi * lphi;

        return done(
            ((((m11 * p2 + m9) * p2 + m7) * p2 + m5) * p2 + m3) * p2 * lphi + lphi,
            sign,
            npio2,
            ellipe,
        );
    }

    let mut t = lphi.tan();
    let mut b = a.sqrt();

    if t.abs() > num!(10.0) {
        let e = one!() / (b * t);
        if e.abs() < num!(10.0) {
            let e = e.atan();
            return done(
                ellipe + m * lphi.sin() * e.sin() - ellipeinc(e, m)?,
                sign,
                npio2,
                ellipe,
            );
        }
    }

    let mut c = m.sqrt();
    let mut a = one!();
    let mut d = one!();
    let mut ee = zero!();
    let mut mod_phi = 0;

    while (c / a).abs() > epsilon!() {
        let temp = b / a;
        lphi = lphi + ((t * temp).atan() + num!(mod_phi) * pi!());
        let denom = one!() - temp * t * t;

        if denom.abs() > num!(10.0) * epsilon!() {
            t = t * (one!() + temp) / denom;
            mod_phi = ((lphi + pi_2!()) / pi!()).to_i32().unwrap();
        } else {
            t = lphi.tan();
            mod_phi = ((lphi - t.atan()) / pi!()).floor().to_i32().unwrap();
        }

        c = (a - b) / two!();
        let temp = (a * b).sqrt();
        a = (a + b) / two!();
        b = temp;
        d = d + d;
        ee = ee + c * lphi.sin();
    }

    done(
        ellipe / ellipk(m)? * (t.atan() + num!(mod_phi) * pi!()) / (d * a) + ee,
        sign,
        npio2,
        ellipe,
    )
}

/// Compute elliptic integral of the second kind for m<0.
#[inline]
fn ellipeinc_neg_m<T: Float>(phi: T, m: T) -> T {
    let mpp = m * phi * phi;

    if -mpp < num!(1e-6) && phi < -m {
        return phi + (mpp * phi * phi / num!(30.0) - mpp * mpp / num!(40.0) - mpp / six!()) * phi;
    }

    if -mpp > num!(1e6) {
        let sm = (-m).sqrt();
        let sp = phi.sin();
        let cp = phi.cos();

        let a = one!() - cp;
        let b1 = (four!() * sp * sm / (one!() + cp)).ln();
        let b = -(half!() + b1) / two!() / m;
        let c = (num!(0.75) + cp / sp / sp - b1) / num!(16.0) / m / m;
        return (a + b + c) * sm;
    }

    let scalef: T;
    let scaled: T;
    let x: T;
    let y: T;
    let z: T;
    if phi > num!(1e-153) && m > num!(-1e200) {
        let s = phi.sin();
        let csc2 = one!() / (s * s);
        scalef = one!();
        scaled = m / three!();
        let phi_tan = phi.tan();
        x = one!() / (phi_tan * phi_tan);
        y = csc2 - m;
        z = csc2;
    } else {
        scalef = phi;
        scaled = mpp * phi / three!();
        x = one!();
        y = one!() - mpp;
        z = one!();
    }

    if x == y && x == z {
        return (scalef + scaled / x) / x.sqrt();
    }

    let a0f = (x + y + z) / three!();
    let a0d = (x + y + three!() * z) / five!();
    let mut af = a0f;
    let mut ad = a0d;

    let mut seriesd = zero!();
    let mut seriesn = one!();

    let mut q = num!(400.0) * (a0f - x).abs().max((a0f - y).abs().max((a0f - z).abs()));

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
        x1 = (x1 + lam) / four!();
        y1 = (y1 + lam) / four!();
        z1 = (z1 + lam) / four!();
        af = (x1 + y1 + z1) / three!();
        ad = (ad + lam) / four!();
        n += 1;
        q = q / four!();
        seriesn = seriesn / four!();
    }

    let two_to_2n = num!(1 << (2 * n));
    let xf = (a0f - x) / af / two_to_2n;
    let yf = (a0f - y) / af / two_to_2n;
    let zf = -(xf + yf);

    let e2f = xf * yf - zf * zf;
    let e3f = xf * yf * zf;

    let mut ret = scalef
        * (one!() - e2f / num!(10.0) + e3f / num!(14.0) + e2f * e2f / num!(24.0)
            - three!() * e2f * e3f / num!(44.0))
        / af.sqrt();

    let xd = (a0d - x) / ad / two_to_2n;
    let yd = (a0d - y) / ad / two_to_2n;
    let zd = -(xd + yd) / three!();

    let e2d = xd * yd - six!() * zd * zd;
    let e3d = (three!() * xd * yd - eight!() * zd * zd) * zd;
    let e4d = three!() * (xd * yd - zd * zd) * zd * zd;
    let e5d = xd * yd * zd * zd * zd;

    ret = ret
        - scaled
            * (one!() - three!() * e2d / num!(14.0)
                + e3d / six!()
                + nine!() * e2d * e2d / num!(88.0)
                - three!() * e4d / num!(22.0)
                - nine!() * e2d * e3d / num!(52.0)
                + three!() * e5d / num!(26.0))
            / two_to_2n
            / ad
            / ad.sqrt();

    ret - three!() * scaled * seriesd
}

#[cfg(test)]
mod tests {
    use crate::compare_test_data;

    use super::*;

    fn ellipeinc_k<T: Float>(inp: &[T]) -> T {
        ellipeinc(inp[0], inp[1] * inp[1]).unwrap()
    }

    #[test]
    fn test_ellipeinc() {
        compare_test_data!(
            "./tests/data/boost/ellipeinc_data.txt",
            ellipeinc_k::<f64>,
            1.1e-15
        );
    }
}
