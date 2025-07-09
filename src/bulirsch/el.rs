/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use num_traits::Float;

use crate::{ellipeinc, ellipf, StrErr};

use super::{cel1, cel2, BulirschConst};

/// Computes [incomplete elliptic integral of the first kind in Bulirsch's form](https://dlmf.nist.gov/19.2.E11_5).
/// ```text
///                 arctan(x)                                                   
///                ⌠                             
///                |               dϑ
/// el1(x, kc)  =  ⎮  ────────────────────────────
///                ⎮      ______________________
///                ⌡   ╲╱ cos²(ϑ) + kc² sin²(ϑ)    
///               0                                                   
/// ```
///
/// ## Parameters
/// - x: tangent of amplitude angle. x ∈ ℝ.
/// - kc: complementary modulus. kc ∈ ℝ, kc ≠ 0.
///
/// The precision of the function can be adjusted by overwriting the trait [super::BulirschConst].
/// The default is set according to the original literature by [Bulirsch](https://doi.org/10.1007/BF02165405) for [f64] and [f32].
///
/// ## Domain
/// - Returns error if kc = 0.
///
/// ## Graph
/// ![Bulirsch's Incomplete Elliptic Integral of the First Kind](https://github.com/p-sira/ellip/blob/main/figures/el1_plot.svg?raw=true)
///
/// [Interactive Plot](https://github.com/p-sira/ellip/blob/main/figures/el1_plot.html)
///
/// # Related Functions
/// With x = tan φ and kc² = 1 - m,
/// - [ellipf](crate::ellipf)(φ, m) = [el1](crate::el1)(x, kc) = [el2](crate::el2)(x, kc, 1, 1)
/// - [el1](crate::el1)(∞, kc) = [cel1](crate::cel1)(kc)
///
/// # Examples
/// ```
/// use ellip::{el1, util::assert_close};
/// use std::f64::consts::FRAC_PI_4;
///
/// assert_close(el1(FRAC_PI_4.tan(), 0.5).unwrap(), 0.8512237490711854, 1e-15);
/// ```
///
/// # References
/// - Bulirsch, Roland. “Numerical Calculation of Elliptic Integrals and Elliptic Functions.” Numerische Mathematik 7, no. 1 (February 1, 1965): 78–90. <https://doi.org/10.1007/BF01397975>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
///
pub fn el1<T: Float + BulirschConst>(x: T, kc: T) -> Result<T, StrErr> {
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
/// Computes [incomplete elliptic integral of the second kind in Bulirsch's form](https://dlmf.nist.gov/19.2.E12).
/// ```text
///                       arctan(x)                                                   
///                      ⌠                             
///                      |              a + b tan²(ϑ)
/// el2(x, kc, a, b)  =  ⎮  ──────────────────────────────────── ⋅ dϑ
///                      ⎮     ________________________________
///                      ⌡  ╲╱ (1 + tan²(ϑ)) (1 + kc² tan²(ϑ))    
///                     0                                                   
/// ```
///
/// ## Parameters
/// - x: tangent of amplitude angle. x ∈ ℝ.
/// - kc: complementary modulus. kc ∈ ℝ, kc ≠ 0.
/// - a ∈ ℝ
/// - b ∈ ℝ
///
/// The precision of the function can be adjusted by overwriting the trait [super::BulirschConst].
/// The default is set according to the original literature by [Bulirsch](https://doi.org/10.1007/BF02165405) for [f64] and [f32].
///
/// ## Domain
/// - Returns error if kc = 0.
///
/// ## Graph
/// ![Bulirsch's Incomplete Elliptic Integral of the Second Kind](https://github.com/p-sira/ellip/blob/main/figures/el2_plot.svg?raw=true)
///
/// [Interactive Plot](https://github.com/p-sira/ellip/blob/main/figures/el2_plot.html)
///
/// # Related Functions
/// With x = tan φ and kc² = 1 - m,
/// - [ellipf](crate::ellipf)(φ, m) = [el1](crate::el1)(x, kc) = [el2](crate::el2)(x, kc, 1, 1)
/// - [ellipeinc](crate::ellipeinc)(φ, m) = [el2](crate::el2)(x, kc, 1, kc²)
/// - [el2](crate::el2)(∞, kc, a, b) = [cel1](crate::cel2)(kc, a, b)
///
/// # Examples
/// ```
/// use ellip::{el2, util::assert_close};
/// use std::f64::consts::FRAC_PI_4;
///
/// assert_close(el2(FRAC_PI_4.tan(), 0.5, 1.0, 1.0).unwrap(), 0.8512237490711854, 1e-15);
/// ```
///
/// # References
/// - Bulirsch, Roland. “Numerical Calculation of Elliptic Integrals and Elliptic Functions.” Numerische Mathematik 7, no. 1 (February 1, 1965): 78–90. <https://doi.org/10.1007/BF01397975>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
///
pub fn el2<T: Float + BulirschConst>(x: T, kc: T, a: T, b: T) -> Result<T, StrErr> {
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

/// Computes [incomplete elliptic integral of the third kind in Bulirsch's form](https://dlmf.nist.gov/19.2.E16).
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
/// ## Parameters
/// - x: tangent of amplitude angle. x ∈ ℝ.
/// - kc: complementary modulus. kc ∈ ℝ, kc ≠ 0.
/// - p ∈ ℝ
///
/// The precision of the function can be adjusted by overwriting the trait [super::BulirschConst].
/// The default is set according to the original literature by [Bulirsch](https://doi.org/10.1007/BF02165405) for [f64] and [f32].
///
/// ## Domain
/// - Returns error if:
///   - kc = 0,
///   - 1 + px² = 0,
///   - or |kc| > 1 for p < 0.
/// - Returns the Cauchy principal value when 1 + px² < 0
///
/// ## Graph
/// ![Bulirsch's Incomplete Elliptic Integral of the Third Kind](https://github.com/p-sira/ellip/blob/main/figures/el3_plot.svg?raw=true)
///
/// [Interactive Plot](https://github.com/p-sira/ellip/blob/main/figures/el3_plot.html)
///
/// # Related Functions
/// With x = tan φ, p = 1 - n, and kc² = 1 - m,
/// - [ellippiinc](crate::ellippiinc)(φ, n, m) = [el3](crate::el3)(x, kc, p)
///
/// # Examples
/// ```
/// use ellip::{el3, util::assert_close};
/// use std::f64::consts::FRAC_PI_4;
///
/// assert_close(el3(FRAC_PI_4.tan(), 0.5, 1.0).unwrap(), 0.8512237490711854, 1e-15);
/// ```
///
/// # References
/// - Bulirsch, R. “Numerical Calculation of Elliptic Integrals and Elliptic Functions. III.” Numerische Mathematik 13, no. 4 (August 1, 1969): 305–15. <https://doi.org/10.1007/BF02165405>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
///
pub fn el3<T: Float + BulirschConst>(x: T, kc: T, p: T) -> Result<T, StrErr> {
    if kc == zero!() {
        return Err("el3: kc must not be zero.");
    }

    if x == zero!() {
        return Ok(zero!());
    }

    // Handle special cases
    let phi = x.atan();
    let m = one!() - kc * kc;
    let n = one!() - p;

    if kc == one!() {
        // A&S 17.7.20:
        if n < one!() {
            let vcr = p.sqrt();
            return Ok((vcr * phi.tan()).atan() / vcr);
        } else {
            // v > 1:
            let vcr = (-p).sqrt();
            let arg = vcr * phi.tan();
            return Ok((arg.ln_1p() - (-arg).ln_1p()) / (two!() * vcr));
        }
    }

    // Didn't improve the accuracy
    // if n == zero!() {
    //     // A&S 17.7.18 & 19
    //     return if m == zero!() {
    //         Ok(phi)
    //     } else {
    //         ellipf(phi, m)
    //     };
    // }

    // This cutpoint is empirical
    if p.abs() <= num!(1e-10) {
        if m == zero!() {
            return Ok(phi.tan());
        }

        // http://functions.wolfram.com/08.06.03.0008.01
        let sp2 = phi.sin() * phi.sin();
        let mut result = (one!() - m * sp2).sqrt() * phi.tan() - ellipeinc(phi, m)?;
        result = result / (one!() - m);
        result = result + ellipf(phi, m)?;
        return Ok(result);
    }

    // // https://dlmf.nist.gov/19.6.E11
    // if phi == zero!() {
    //     return Ok(zero!());
    // }

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
        for k in 0..nd - 2 {
            rr[k] = s;
            pm = pm * t * ra[k];
            s = s * p + pm;
        }
        s = s * zd;
        u = s;
        bo = false;
        for k in (0..nd - 2).rev() {
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
        return Err("el3: 1 + px² cannot be zero.");
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
            k = if r < zero!() { -1 } else { 1 };
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
            m += if h < zero!() { -1 } else { 1 };
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
            let m = 1.0 - kc * kc;
            let phi = x.atan();
            let f = ellipf(phi, m).unwrap();

            assert_close!(f, el1(x, kc).unwrap(), 1.2e-14);
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
            let m = 1.0 - kc * kc;
            let phi = x.atan();
            let f = ellipf(phi, m).unwrap();
            let e = ellipeinc(phi, m).unwrap();
            let d = ellipdinc(phi, m).unwrap();

            assert_close!(f, el2(x, kc, 1.0, 1.0).unwrap(), 1.2e-14);
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
            let m = 1.0 - kc * kc;
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

        // Test computed values from the reference
        // Bulirsch, “Numerical Calculation of Elliptic Integrals and Elliptic Functions III.”
        fn test_reference(x: f64, kc: f64, p: f64, expected: f64) {
            assert_close!(expected, el3(x, kc, p).unwrap(), 2.0 * f64::EPSILON);
        }

        test_reference(1.3, 0.11, 4.21, 6.6220785847015254e-1);
        test_reference(1.3, 0.11, 0.82, 1.1307046442074609);
        test_reference(1.3, 0.92, 0.71, 1.0058286266977115);
        test_reference(1.3, 0.92, 0.23, 1.1884070823345123);
        test_reference(1.3, 0.12, -0.11, 1.7259650355348878);
        test_reference(1.3, 0.12, -2.11, 2.4416814520721179e-1);
        test_reference(1.3, 0.40, 0.1600001, 1.4004165258366944);
        test_reference(1.3, 1.0e-10, 0.82, 1.1341505395282723);
        test_reference(1.3e-10, 1.0e-10, 1.0e-10, 1.3e-10);
        test_reference(1.6, 1.90, 9.81, 3.8572324379967252e-1);
        test_reference(1.6, 1.90, 1.22, 7.6656179311956402e-1);
        test_reference(1.6, 1.90, 0.87, 8.3210591112618096e-1);
        test_reference(1.6, 1.90, 0.21, 1.0521272221906806);
        test_reference(1.6, 1.90, -0.21, 1.4730439889554361);
        test_reference(1.6, 1.90, -4.30, 2.5467519341311686e-1);
        test_reference(1.6, 1.01e1, -1.0e-5, 3.9501709882649139e-1);
        test_reference(1.6, 1.50, 2.24999, 7.0057431688357934e-1);
        test_reference(1.6, 1e10, 1.20, 2.3734774669772208e-9);
        test_reference(-1.6, 1e10, 1.20, -2.3734774669772208e-9);
        test_reference(1.0, 0.31, 9.90e-2, 1.0903577921777398);
    }
}
