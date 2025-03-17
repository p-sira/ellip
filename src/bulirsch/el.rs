/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use num_traits::Float;

use super::{cel1, cel2, BulirschConst};

// Reference: Bulirsch, “Numerical Calculation of Elliptic Integrals and Elliptic Functions.”
/// Compute [incomplete elliptic integral of the first kind in Bulirsch's form](https://dlmf.nist.gov/19.2.E11_5).
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
/// Note that x = tan φ and kc² = mc = 1 - m. The precision can be adjusted by overwriting the trait [super::BulirschConst].
///
/// # Examples
/// ```
/// use ellip::{el1, util::assert_close};
/// use std::f64::consts::FRAC_PI_4;
///
/// assert_close(el1(FRAC_PI_4.tan(), 0.5).unwrap(), 0.8512237490711854, 1e-15);
/// ```
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
/// Compute [incomplete elliptic integral of the second kind in Bulirsch's form](https://dlmf.nist.gov/19.2.E12).
/// ```text
///                       arctan(x)                                                   
///                      ⌠                             
///                      |              a + b tan²(ϑ)
/// el2(x, kc, a, b)  =  ⎮  ──────────────────────────────────── ⋅ dϑ
///                      ⎮     ________________________________
///                      ⌡  ╲╱ (1 + tan²(ϑ)) (1 + kc² tan²(ϑ))    
///                     0                                                   
/// where kc ≠ 0
/// ```
///
/// Note that x = tan φ and kc² = mc = 1 - m. The precision can be adjusted by overwriting the trait [super::BulirschConst].
///
/// # Examples
/// ```
/// use ellip::{el2, util::assert_close};
/// use std::f64::consts::FRAC_PI_4;
///
/// assert_close(el2(FRAC_PI_4.tan(), 0.5, 1.0, 1.0).unwrap(), 0.8512237490711854, 1e-15);
/// ```
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

// Reference: Bulirsch, “Numerical Calculation of Elliptic Integrals and Elliptic Functions. III”
/// Compute [incomplete elliptic integral of the third kind in Bulirsch's form](https://dlmf.nist.gov/19.2.E16).
/// ```text
///                    arctan(x)                                                   
///                   ⌠                             
///                   |                          dϑ
/// el3(x, kc, p)  =  ⎮  ───────────────────────────────────────────────────
///                   ⎮                              ______________________
///                   ⌡  ( cos²(ϑ) + p sin²(ϑ) ) ⋅ ╲╱ cos²(ϑ) + kc² sin²(ϑ)   
///                  0                                                   
/// ```
///
/// Note that x = tan φ and kc² = mc = 1 - m. The precision can be adjusted by overwriting the trait [super::BulirschConst].
///
/// # Examples
/// ```
/// use ellip::{el3, util::assert_close};
/// use std::f64::consts::FRAC_PI_4;
///
/// assert_close(el3(FRAC_PI_4.tan(), 0.5, 1.0).unwrap(), 0.8512237490711854, 1e-15);
/// ```
pub fn el3<T: Float + BulirschConst>(x: T, kc: T, p: T) -> Result<T, &'static str> {
    if x == zero!() {
        return Ok(zero!());
    }

    // real
    let mut c;
    let mut d;
    let mut de;
    let mut e;
    let mut f;
    let mut fa;
    let mut g;
    let mut h;
    let mut hh;

    let mut pm;
    let mut pz;
    let mut q;
    let mut r;
    let mut s;
    let mut t;
    let mut u;
    let mut v;
    let mut w;
    let mut y;
    let mut ye = zero!();
    let mut z;
    let zd;

    // int
    let mut k = 0;
    let mut l;
    let mut m;
    let mut n;
    let nd = T::D as usize - 2;

    // bool
    let mut bo;
    let mut bk = false;

    hh = x * x;
    f = p * hh;
    s = if kc == zero!() {
        T::ca() / (one!() + x.abs())
    } else {
        kc
    };
    t = s * s;
    pm = t * half!();
    e = hh * t;
    z = f.abs();
    r = p.abs();
    h = one!() + hh;

    // small
    if e < tenth!() && z < tenth!() && t < one!() && r < one!() {
        let (rb, ra): (Vec<T>, Vec<T>) = (2..=nd)
            .map(|k| {
                let rb_k = half!() / num!(k);
                (rb_k, one!() - rb_k)
            })
            .unzip();

        zd = half!() / (num!(nd) + one!());
        s = p + pm;

        let mut rr: Vec<T> = vec![zero!(); nd - 2];
        for k in 0..=nd - 2 {
            rr[k] = s;
            pm = pm * t * ra[k];
            s = s * p + pm;
        }
        s = s * zd;
        u = s;
        bo = false;
        for k in (0..=nd - 2).rev() {
            u = u + (rr[k] - u) * rb[k];
            bo = !bo;
            let v = if bo { -u } else { u };
            s = s * hh + v;
        }

        if bo {
            s = -s;
        }
        u = half!() * (u + one!());
        return Ok((u - s * h) * h.sqrt() * x + u * x.asinh());
    }

    w = one!() + f;
    if w == zero!() {
        return Err("fail");
    }

    let p1 = if p == zero!() { T::cb() / hh } else { p };
    s = s.abs();
    y = x.abs();
    g = p1 - one!();
    if g == zero!() {
        g = T::cb();
    }
    f = p1 - t;
    if f == zero!() {
        f = T::cb() * t;
    }
    let am = one!() - t;
    let ap = one!() + e;
    r = p1 * h;
    fa = g / (f * p1);
    bo = fa > zero!();
    fa = fa.abs();
    pz = (g * f).abs();
    de = pz.sqrt();
    q = p1.abs().sqrt();

    if pm > half!() {
        pm = half!()
    } else {
        pm = p1 - pm
    };

    if pm >= zero!() {
        u = (r * ap).sqrt();
        v = y * de;
        if g < zero!() {
            v = -v;
        }
        d = one!() / q;
        c = one!();
    } else {
        u = (h * ap * pz).sqrt();
        ye = y * q;
        v = am * ye;
        q = -de / g;
        d = -am / de;
        c = zero!();
        pz = ap - r;
    }

    if bo {
        r = v / u;
        z = one!();
        k = 1;
        if pm < zero!() {
            h = y * (h / ap / fa).sqrt();
            h = one!() / h - h;
            z = h - r - r;
            r = two!() + r * h;
            if r == zero!() {
                r = T::cb();
            }
            if z == zero!() {
                z = h * T::cb();
            }
            z = r / z;
            r = z;
            w = pz;
        }
        u = u / w;
        v = v / w;
    } else {
        t = u + v.abs();
        bk = true;
        if p1 < zero!() {
            de = v / pz;
            ye = u * ye;
            ye = ye + ye;
            u = t / pz;
            v = (-f - g * e) / t;
            t = pz * w.abs();
            z = (hh * r * f - g * ap + ye) / t;
            ye = ye / t;
        } else {
            de = v / w;
            ye = zero!();
            u = (e + p1) / t;
            v = t / w;
            z = one!();
        }
        if s > one!() {
            h = u;
            u = v;
            v = h;
        }
    }

    y = one!() / y;
    e = s;
    n = 1;
    t = one!();
    l = 0;
    m = 0;

    loop {
        y = y - e / y;
        if y == zero!() {
            y = e.sqrt() * T::cb();
        }
        f = c;
        c = d / q + c;
        g = e / q;
        d = f * g + d;
        d = d + d;
        q = g + q;
        g = t;
        t = s + t;
        n = n + n;
        m = m + m;
        if bo {
            if z < zero!() {
                m += k;
            }
            k = r.signum().to_i32().expect("el3: Cannot extract sign of r.");
            h = e / (u * u + v * v);
            u = u * (one!() + h);
            v = v * (one!() - h);
        } else {
            r = u / v;
            h = z * r;
            z = h * z;
            hh = e / v;

            de = de / u;
            ye = ye * (h + one!() / h) + de * (one!() + r);
            de = de * (u - hh);
            bk = ye.abs() < one!();

            // Removed crack function
        }
        if (g - s).abs() > T::ca() * g {
            if bo {
                g = (one!() / r - r) * half!();
                hh = u + v * g;
                h = g * u - v;
                if hh == zero!() {
                    hh = u * T::cb();
                }
                if h == zero!() {
                    h = v * T::cb();
                }
                z = r * h;
                r = hh / h;
            } else {
                u = u + e / u;
                v = v + hh;
            }
            s = e.sqrt();
            s = s + s;
            e = s * t;
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
    e = (t / y).atan() + pi!() * num!(l);
    e = e * (c * t + d) / (t * (t + q));

    if bo {
        h = v / (t + u);
        z = one!() - r * h;
        h = r + h;
        if z == zero!() {
            z = T::cb();
        }
        if z < zero!() {
            m += h.signum().to_i32().expect("el3: Cannot extract sign of h.");
        }
        s = (h / z).atan() + num!(m) * pi!();
    } else {
        s = if bk {
            ye.asinh()
        } else {
            z.ln() + num!(m) * ln_2!()
        };
        s = s * half!();
    }
    e = (e + fa.sqrt() * s) / num!(n);

    Ok(if x > zero!() { e } else { -e })
}

#[cfg(test)]
mod tests {
    use itertools::iproduct;

    use super::*;
    use crate::{assert_close, ellipdinc, ellipeinc, ellipf, ellippiinc, test_util::linspace};

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
            assert_close!(e, el2(x, kc, 1.0, kc * kc).unwrap(), 7e-14);
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

    #[test]
    fn test_el3() {
        fn _test(x: f64, kc: f64, p: f64) {
            let k = (1.0 - kc * kc).sqrt();
            let m = k * k;
            let phi = x.atan();
            let n = 1.0 - p;
            let ellippi = ellippiinc(phi, n, m).unwrap();

            assert_close!(ellippi, el3(x, kc, p).unwrap(), 4.6e-15);
        }

        let linsp_x = linspace(10.0, 10.0, 15);
        let linsp_kc = [
            linspace(-1.0 + 1e-3, -1e-3, 50),
            linspace(1e-3, 1.0 - 1e-3, 50),
        ]
        .concat();
        let linsp_p = linspace(1e-3, 1.0 - 1e-3, 10);

        iproduct!(linsp_x, linsp_kc, linsp_p).for_each(|(x, kc, p)| _test(x, kc, p));
    }
}
