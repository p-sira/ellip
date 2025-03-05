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
/// where 0 ≤ m sin²φ ≤ 1
/// ```
///
/// Note that some mathematical references use the parameter k for the function,
/// where k² = m.
///
/// # Examples
/// ```
/// use ellip::{ellipf, util::assert_close};
/// use std::f64::consts::FRAC_PI_4;
///
/// assert_close(ellipf(FRAC_PI_4, 0.5).unwrap(), 0.826017876249245, 1e-15);
/// ```
pub fn ellipf<T: Float>(phi: T, m: T) -> Result<T, &'static str> {
    let ms2p =  m * (phi.sin() * phi.sin());
    if ms2p > one!() || ms2p < zero!() {
        return Err("ellipf: m sin²φ must satisfy: 0 ≤ m sin²φ ≤ 1.");
    }

    if phi.is_infinite() || m.is_infinite() {
        if phi.is_infinite() {
            return Ok(inf!());
        }

        if m.is_infinite() {
            return Ok(zero!());
        }

        return Err("ellipf: m or φ must be finite");
    }

    if m == zero!() {
        return Ok(phi);
    }

    if m == one!() {
        if phi.abs() >= pi_2!() {
            return Ok(inf!());
        }
        return Ok(phi.tan().sinh().atanh());
    }

    let mut npio2 = (phi / pi_2!()).floor();
    if (npio2.abs() % two!()) == one!() {
        npio2 = npio2 + one!();
    }

    let mut phi = phi;
    let mut k;
    if npio2 != zero!() {
        k = ellipk(m)?;
        phi = phi - npio2 * pi_2!();
    } else {
        k = zero!();
    }

    let sign = if phi < zero!() {
        phi = -phi;
        -one!()
    } else {
        one!()
    };

    fn done<T: Float>(val: T, sign: T, npio2: T, k: T) -> Result<T, &'static str> {
        Ok(sign * val + npio2 * k)
    }

    let a = one!() - m;
    if a > one!() {
        return done(ellipf_neg_m(phi, m), sign, npio2, k);
    }

    let mut b = a.sqrt();
    let mut t = phi.tan();
    if t.abs() > num!(10.0) {
        let mut e = one!() / (b * t);
        if e.abs() < num!(10.0) {
            e = e.atan();
            if npio2 == zero!() {
                k = ellipk(m)?;
            }
            return done(k - ellipf(e, m)?, sign, npio2, k);
        }
    }

    let mut a = one!();
    let mut c = m.sqrt();
    let mut d = one!();
    let mut mod_phi = zero!();

    while (c / a).abs() > epsilon!() {
        let temp = b / a;
        phi = phi + (t * temp).atan() + mod_phi * pi!();
        let denom = one!() - temp * t * t;
        if denom.abs() > num!(10.0) * epsilon!() {
            t = t * (one!() + temp) / denom;
            mod_phi = num!(((phi + pi_2!()) / pi!()).to_i32().unwrap());
        } else {
            t = phi.tan();
            mod_phi = num!(((phi - t.atan()) / pi!()).floor().to_i32().unwrap());
        }
        c = (a - b) / two!();
        let temp = (a * b).sqrt();
        a = (a + b) / two!();
        b = temp;
        d = d * two!();
    }

    done((t.atan() + mod_phi * pi!()) / (d * a), sign, npio2, k)
}

/// Compute elliptic integral of the first kind for m<0.
#[inline]
fn ellipf_neg_m<T: Float>(phi: T, m: T) -> T {
    let mpp = m * phi * phi;

    if -mpp < num!(1e-6) && phi < -m {
        return phi
            + (-mpp * phi * phi / num!(30.0) + three!() * mpp * mpp / num!(40.0) + mpp / six!())
                * phi;
    }

    if -mpp > num!(4e7) {
        let sm = (-m).sqrt();
        let sp = phi.sin();
        let cp = phi.cos();

        let a = (four!() * sp * sm / (one!() + cp)).ln();
        let b = -(one!() + cp / (sp * sp) - a) / (four!() * m);
        return (a + b) / sm;
    }

    let (scale, x, y, z) = if phi > num!(1e-153) && m > num!(-1e305) {
        let s = phi.sin();
        let phi_tan = phi.tan();
        let csc2 = one!() / (s * s);
        (one!(), one!() / (phi_tan * phi_tan), csc2 - m, csc2)
    } else {
        (phi, one!(), one!() - m * phi * phi, one!())
    };

    if x == y && x == z {
        return scale / x.sqrt();
    }

    let a0 = (x + y + z) / three!();
    let mut a = a0;
    let mut x1 = x;
    let mut y1 = y;
    let mut z1 = z;
    let mut q = num!(400.0) * (a0 - x).abs().max((a0 - y).abs()).max((a0 - z).abs());
    let mut n = 0;

    while q > a.abs() && n <= 100 {
        let sx = x1.sqrt();
        let sy = y1.sqrt();
        let sz = z1.sqrt();
        let lam = sx * sy + sx * sz + sy * sz;
        x1 = (x1 + lam) / four!();
        y1 = (y1 + lam) / four!();
        z1 = (z1 + lam) / four!();
        a = (x1 + y1 + z1) / three!();
        n += 1;
        q = q / four!();
    }

    let two_to_2n = num!(1u32 << (2 * n));
    let x = (a0 - x) / a / two_to_2n;
    let y = (a0 - y) / a / two_to_2n;
    let z = -(x + y);

    let e2 = x * y - z * z;
    let e3 = x * y * z;

    scale
        * (one!() - e2 / num!(10.0) + e3 / num!(14.0) + e2 * e2 / num!(24.0)
            - three!() * e2 * e3 / num!(44.0))
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
        compare_test_data!("./tests/data/boost/ellipf_data.txt", ellipf_k, 4.9e-16);
    }
}
