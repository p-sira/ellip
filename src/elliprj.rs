/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use std::mem::swap;

use crate::{elliprc, elliprd, elliprf};

pub fn elliprj(x: f64, y: f64, z: f64, p: f64) -> Result<f64, &'static str> {
    if x.min(y).min(z) < 0.0 || (y + z).min(x + y).min(x + z) == 0.0 {
        return Err("elliprj: x, y, and z must be non-negative, and at most one can be zero.");
    }
    if p == 0.0 {
        return Err("elliprj: p must be non-zero");
    }

    // for p < 0, the integral is singular, return Cauchy principal value
    if p < 0.0 {
        // We must ensure that x < y < z.
        // Since the integral is symmetrical in x, y and z
        // we can just permute the values:
        let mut x = x;
        let mut y = y;
        let mut z = z;

        if x > y {
            swap(&mut x, &mut y);
        }
        if y > z {
            swap(&mut y, &mut z);
        }
        if x > y {
            swap(&mut x, &mut y);
        }

        let q = -p;
        let p = (z * (x + y + q) - x * y) / (z + q);

        let value = (p - z) * _elliprj(x, y, z, p)? - 3.0 * elliprf(x, y, z)?
            + 3.0 * ((x * y * z) / (x * y + p * q)).sqrt() * elliprc(x * y + p * q, p * q)?;
        Ok(value / (z + q))
    } else {
        _elliprj(x, y, z, p)
    }
}

/// Calculate RC(1, 1 + x)
fn elliprc1p(y: f64) -> Result<f64, &'static str> {
    // We can skip this check since the call from _elliprj already did the check.
    // if y == -1.0 {
    //     return Err("elliprc1p: y cannot be -1.0.");
    // }

    // for 1 + y < 0, the integral is singular, return Cauchy principal value
    if y < -1.0 {
        Ok((1.0 / -y).sqrt() * elliprc(-y, -1.0 - y)?)
    } else if y == 0.0 {
        Ok(1.0)
    } else if y > 0.0 {
        Ok(y.sqrt().atan() / y.sqrt())
    } else if y > -0.5 {
        let arg = (-y).sqrt();
        Ok(((1.0 + arg).ln() - (1.0 - arg).ln()) / (2.0 * (-y).sqrt()))
    } else {
        Ok(((1.0 + (-y).sqrt()) / (1.0 + y).sqrt()).ln() / (-y).sqrt())
    }
}

fn _elliprj(mut x: f64, y: f64, mut z: f64, p: f64) -> Result<f64, &'static str> {
    // Special cases
    // https://dlmf.nist.gov/19.20#iii
    if x == y {
        if x == z {
            if x == p {
                // RJ(x,x,x,x)
                return Ok(1.0 / (x * x.sqrt()));
            } else {
                // RJ(x,x,x,p)
                return Ok((3.0 / (x - p)) * (elliprc(x, p)? - 1.0 / x.sqrt()));
            }
        } else {
            // RJ(x,x,z,p)
            swap(&mut x, &mut z);
            // RJ(x,y,y,p)
            // Fall through to next if block.
        }
    }

    if y == z {
        if y == p {
            // RJ(x,y,y,y)
            return elliprd(x, y, y);
        }
        // This prevents division by zero.
        if p.max(y) / p.min(y) > 1.2 {
            // RJ(x,y,y,p)
            return Ok((3.0 / (p - y)) * (elliprc(x, y)? - elliprc(x, p)?));
        }
    }

    if z == p {
        // RJ(x,y,z,z)
        return elliprd(x, y, z);
    }

    let mut xn = x;
    let mut yn = y;
    let mut zn = z;
    let mut pn = p;
    let mut an = (x + y + z + 2.0 * p) / 5.0;
    let a0 = an;
    let mut delta = (p - x) * (p - y) * (p - z);
    let q = (f64::EPSILON / 5.0).powf(-1.0 / 8.0)
        * (an - x)
            .abs()
            .max((an - y).abs())
            .max((an - z).abs())
            .max((an - p).abs());

    let mut fmn = 1.0; // 4^-n
    let mut rc_sum = 0.0;

    for _ in 0..100 {
        // max iterations
        let rx = xn.sqrt();
        let ry = yn.sqrt();
        let rz = zn.sqrt();
        let rp = pn.sqrt();
        let dn = (rp + rx) * (rp + ry) * (rp + rz);
        let en = delta / (dn * dn);

        if en < -0.5 && en > -1.5 {
            let b = 2.0 * rp * (pn + rx * (ry + rz) + ry * rz) / dn;
            rc_sum += fmn / dn * elliprc(1.0, b)?;
        } else {
            rc_sum += fmn / dn * elliprc1p(en)?;
        }

        let lambda = rx * ry + rx * rz + ry * rz;
        an = (an + lambda) / 4.0;
        fmn /= 4.0;

        if fmn * q < an {
            break;
        }

        xn = (xn + lambda) / 4.0;
        yn = (yn + lambda) / 4.0;
        zn = (zn + lambda) / 4.0;
        pn = (pn + lambda) / 4.0;
        delta /= 64.0;
    }

    let x = fmn * (a0 - x) / an;
    let y = fmn * (a0 - y) / an;
    let z = fmn * (a0 - z) / an;
    let p = (-x - y - z) / 2.0;
    let xyz = x * y * z;
    let p2 = p * p;
    let p3 = p2 * p;

    let e2 = x * y + x * z + y * z - 3.0 * p2;
    let e3 = xyz + 2.0 * e2 * p + 4.0 * p3;
    let e4 = (2.0 * xyz + e2 * p + 3.0 * p3) * p;
    let e5 = xyz * p2;

    let result = fmn
        * an.powf(-1.5)
        * (1.0 - 3.0 * e2 / 14.0 + e3 / 6.0 + 9.0 * e2 * e2 / 88.0
            - 3.0 * e4 / 22.0
            - 9.0 * e2 * e3 / 52.0
            + 3.0 * e5 / 26.0
            - e2 * e2 * e2 / 16.0
            + 3.0 * e3 * e3 / 40.0
            + 3.0 * e2 * e4 / 20.0
            + 45.0 * e2 * e2 * e3 / 272.0
            - 9.0 * (e3 * e4 + e2 * e5) / 68.0);

    Ok(result + 6.0 * rc_sum)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::compare_test_data;

    fn elliprj_wrapper(inp: &[f64]) -> f64 {
        elliprj(inp[0], inp[1], inp[2], inp[3]).unwrap()
    }

    #[test]
    fn test_elliprj() {
        compare_test_data!(
            "./tests/data/boost/ellint_rj_data.txt",
            elliprj_wrapper,
            3e-8, // Relatively low precision
            5e-25
        );
    }
}
