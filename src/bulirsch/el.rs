/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use num_traits::Float;

use crate::{
    bulirsch::constants::{BulirschConst, DefaultPrecision},
    cel1, cel2,
    crate_util::{case, check, declare, let_mut},
    ellipeinc, ellipf, ellippi, ellippiinc, StrErr,
};

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
/// ## Domain
/// - Returns error if kc = 0.
///
/// ## Graph
/// ![Bulirsch's Incomplete Elliptic Integral of the First Kind](https://github.com/p-sira/ellip/blob/main/figures/el1_plot.svg?raw=true)
///
/// [Interactive Plot](https://github.com/p-sira/ellip/blob/main/figures/el1_plot.html)
///
/// ## Special Cases
/// - el1(0, kc) = 0
/// - el1(∞, kc) = cel1(kc)
/// - el1(x, ∞) = 0
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
/// # Notes
/// The default precision of the function is set according to the original literature by [Bulirsch](https://doi.org/10.1007/BF02165405)
/// for [f64]. The precision can be modified in the function [_el1] (requires `unstable` feature flag).
///
/// # References
/// - Bulirsch, Roland. “Numerical Calculation of Elliptic Integrals and Elliptic Functions.” Numerische Mathematik 7, no. 1 (February 1, 1965): 78–90. <https://doi.org/10.1007/BF01397975>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
pub fn el1<T: Float>(x: T, kc: T) -> Result<T, StrErr> {
    _el1::<T, DefaultPrecision>(x, kc)
}

/// Computes [el1]. Control the precision using [BulirschConst].
/// <div class="warning">⚠️ Unstable feature. May subject to changes.</div>
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub fn _el1<T: Float, C: BulirschConst<T>>(x: T, kc: T) -> Result<T, StrErr> {
    let ans = el1_unchecked::<T, C>(x, kc);
    if ans.is_finite() {
        return Ok(ans);
    }
    check!(@nan, el1, [x, kc]);
    check!(@zero, el1, [kc]);
    case!(kc == inf!(), T::zero());
    case!(x == T::zero(), T::zero());
    if x == inf!() {
        return cel1(kc);
    }
    Err("el1: Failed to converge.")
}

/// Unsafe version of [el1].
/// <div class="warning">⚠️ Unstable feature. May subject to changes.</div>
/// Undefined behavior with invalid arguments and edge cases.
/// # Known Invalid Cases
/// - kc = 0
/// - x = 0
/// - kc = ∞ or x = ∞
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub fn el1_unchecked<T: Float, C: BulirschConst<T>>(x: T, kc: T) -> T {
    declare!(mut [y = x.recip().abs(), kc = kc.abs(), m = T::one(), l = 0, e, g]);

    for _ in 0..N_MAX_ITERATIONS {
        e = m * kc;
        g = m;
        m = kc + m;
        y = -e / y + y;

        if y == 0.0 {
            y = e.sqrt() * C::cb();
        }

        if (g - kc).abs() > C::ca() * g {
            kc = e.sqrt() * 2.0;
            l *= 2;
            if y < 0.0 {
                l += 1;
            }
            continue;
        }

        if y < 0.0 {
            l += 1;
        }

        return x.signum() * ((m / y).atan() + pi!() * T::from(l).unwrap()) / m;
    }
    nan!()
}

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
/// ## Domain
/// - Returns error if kc = 0.
///
/// ## Graph
/// ![Bulirsch's Incomplete Elliptic Integral of the Second Kind](https://github.com/p-sira/ellip/blob/main/figures/el2_plot.svg?raw=true)
///
/// [Interactive Plot](https://github.com/p-sira/ellip/blob/main/figures/el2_plot.html)
///
/// ## Special Cases
/// - el2(0, kc, a, b) = 0
/// - el2(x, kc, 0, 0) = 0
/// - el2(∞, kc, a, b) = cel2(kc, a, b)
///
/// # Related Functions
/// With x = tan φ and kc² = 1 - m,
/// - [ellipf](crate::ellipf)(φ, m) = [el1](crate::el1)(x, kc) = [el2](crate::el2)(x, kc, 1, 1)
/// - [ellipeinc](crate::ellipeinc)(φ, m) = [el2](crate::el2)(x, kc, 1, kc²)
/// - [el2](crate::el2)(∞, kc, a, b) = [cel2](crate::cel2)(kc, a, b)
///
/// # Examples
/// ```
/// use ellip::{el2, util::assert_close};
/// use std::f64::consts::FRAC_PI_4;
///
/// assert_close(el2(FRAC_PI_4.tan(), 0.5, 1.0, 1.0).unwrap(), 0.8512237490711854, 1e-15);
/// ```
///
/// # Notes
/// The default precision of the function is set according to the original literature by [Bulirsch](https://doi.org/10.1007/BF02165405)
/// for [f64]. The precision can be modified in the function [_el2] (requires `unstable` feature flag).
///
/// # References
/// - Bulirsch, Roland. “Numerical Calculation of Elliptic Integrals and Elliptic Functions.” Numerische Mathematik 7, no. 1 (February 1, 1965): 78–90. <https://doi.org/10.1007/BF01397975>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
pub fn el2<T: Float>(x: T, kc: T, a: T, b: T) -> Result<T, StrErr> {
    _el2::<T, DefaultPrecision>(x, kc, a, b)
}

/// Computes [el2]. Control the precision using [BulirschConst].
/// <div class="warning">⚠️ Unstable feature. May subject to changes.</div>
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub fn _el2<T: Float, C: BulirschConst<T>>(x: T, kc: T, a: T, b: T) -> Result<T, StrErr> {
    let ans = el2_unchecked::<T, C>(x, kc, a, b);
    if ans.is_finite() {
        return Ok(ans);
    }
    check!(@nan, el2, [x, kc, a, b]);
    check!(@zero, el2, [kc]);
    case!(x == T::zero(), T::zero());
    if x == inf!() {
        // phi = π/2
        return cel2(kc, a, b);
    }

    Err("el2: Failed to converge.")
}

/// Unsafe version of [el2].
/// <div class="warning">⚠️ Unstable feature. May subject to changes.</div>
/// Undefined behavior with invalid arguments and edge cases.
/// # Known Invalid Cases
/// - kc = 0
/// - x = 0
/// - x = ∞
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub fn el2_unchecked<T: Float, C: BulirschConst<T>>(x: T, kc: T, a: T, b: T) -> T {
    let_mut!(b);
    declare!(mut [c = x * x, d = T::one() + c, p = ((T::one() + kc * kc * c) / d).sqrt()]);

    d = x / d;
    c = d / (p * 2.0);
    let z = a - b;
    let mut i = a;
    let mut a = (b + a) / 2.0;
    declare!(mut [y = x.recip().abs(), f = T::zero(), l = 0, m = T::one(), kc = kc.abs(), e, g]);

    for _ in 0..N_MAX_ITERATIONS {
        b = i * kc + b;
        e = m * kc;
        g = e / p;
        d = f * g + d;
        f = c;
        i = a;
        p = g + p;
        c = (d / p + c) / 2.0;
        g = m;
        m = kc + m;
        a = (b / m + a) / 2.0;
        y = -e / y + y;

        if y == 0.0 {
            y = e.sqrt() * C::cb();
        }

        if (g - kc).abs() > C::ca() * g {
            kc = e.sqrt() * 2.0;
            l *= 2;
            if y < 0.0 {
                l += 1;
            }
            continue;
        }

        if y < 0.0 {
            l += 1;
        }

        let e = x.signum() * ((m / y).atan() + pi!() * T::from(l).unwrap()) * a / m;
        return e + c * z;
    }
    nan!()
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
/// ## Special Cases
/// - el3(0, kc, p) = 0
/// - el3(∞, kc, p) = cel(kc, p, 1, 1) = Π(1-p, 1-kc²)
///
/// # Related Functions
/// With x = tan φ, kc² = 1 - m, and p = 1 - n,
/// - [ellippiinc](crate::ellippiinc)(φ, n, m) = [el3](crate::el3)(x, kc, p)
/// - [el3](crate::el3)(∞, kc, p) = [cel](crate::cel)(kc, p, 1, 1) = [ellippi](crate::ellippi)(n, m)
///
/// # Examples
/// ```
/// use ellip::{el3, util::assert_close};
/// use std::f64::consts::FRAC_PI_4;
///
/// assert_close(el3(FRAC_PI_4.tan(), 0.5, 1.0).unwrap(), 0.8512237490711854, 1e-15);
/// ```
///
/// # Notes
/// The default precision of the function is set according to the original literature by [Bulirsch](https://doi.org/10.1007/BF02165405)
/// for [f64]. The precision can be modified in the function [_el3] (requires `unstable` feature flag).
///
/// # References
/// - Bulirsch, R. “Numerical Calculation of Elliptic Integrals and Elliptic Functions. III.” Numerische Mathematik 13, no. 4 (August 1, 1969): 305–15. <https://doi.org/10.1007/BF02165405>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
pub fn el3<T: Float>(x: T, kc: T, p: T) -> Result<T, StrErr> {
    _el3::<T, DefaultPrecision>(x, kc, p)
}

/// Computes [el3]. Control the precision using [BulirschConst].
/// <div class="warning">⚠️ Unstable feature. May subject to changes.</div>
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub fn _el3<T: Float, C: BulirschConst<T>>(x: T, kc: T, p: T) -> Result<T, StrErr> {
    let m = 1.0 - kc * kc;
    let n = 1.0 - p;

    if kc.abs() < epsilon!() {
        return Err("el3: kc must not be zero.");
    }

    if x == 0.0 {
        return Ok(0.0);
    }

    // Handle special cases
    if x == inf!() {
        return ellippi(n, m);
    }

    if kc == 1.0 {
        // A&S 17.7.20:
        if n < 1.0 {
            let vcr = p.sqrt();
            return Ok((vcr * x).atan() / vcr);
        } else {
            // v > 1:
            let vcr = (-p).sqrt();
            let arg = vcr * x;
            return Ok((arg.ln_1p() - (-arg).ln_1p()) / (2.0 * vcr));
        }
    }

    // Didn't improve the accuracy
    // if n = 0.0 {
    //     // A&S 17.7.18 & 19
    //     return if m = 0.0 {
    //         Ok(phi)
    //     } else {
    //         ellipf(phi, m)
    //     };
    // }

    // This cutpoint is empirical
    let phi = x.atan();
    if p.abs() <= 1e-10 && m != 1.0 {
        if m == 0.0 {
            return Ok(x);
        }

        // http://functions.wolfram.com/08.06.03.0008.01
        let sp2 = phi.sin() * phi.sin();
        let mut result = (1.0 - m * sp2).sqrt() * x - ellipeinc(phi, m)?;
        result = result / (1.0 - m);
        result = result + ellipf(phi, m)?;
        return Ok(result);
    }

    if kc.abs() < C::lim_kc_p() && p > C::lim_kc_p() {
        return ellippiinc(phi, n, m);
    }

    let x_abs = x.abs();

    // real
    declare!(mut [c, d, de, e, f, fa, g, h, hh]);
    declare!(mut [pm, pz, q, r, s, t, u, v, w, y, ye = T::zero(), z]);
    let zd;

    // int
    declare!(mut [k = 0, l, mm, nn]);

    // bool
    declare!(mut [bo, bk = false]);

    hh = x * x;
    f = p * hh;
    s = if kc == 0.0 {
        C::ca() / (1.0 + x_abs)
    } else {
        kc
    };
    t = s * s;
    pm = t * 0.5;
    e = hh * t;
    z = f.abs();
    r = p.abs();
    h = 1.0 + hh;

    // small
    if e < 0.1 && z < 0.1 && t < 1.0 && r < 1.0 {
        declare!(mut [rb, ra, rr] = [T::zero(); MAX_ND]);

        for k in 2..=C::ND {
            let k_float = T::from(k).unwrap();
            rb[k - 2] = 0.5 / k_float;
            ra[k - 2] = 1.0 - rb[k - 2];
        }

        zd = 0.5 / (T::from(C::ND).unwrap() + 1.0);
        s = p + pm;

        for k in 0..C::ND - 2 {
            rr[k] = s;
            pm = pm * t * ra[k];
            s = s * p + pm;
        }
        s = s * zd;
        u = s;
        bo = false;
        for k in (0..C::ND - 2).rev() {
            u = u + (rr[k] - u) * rb[k];
            bo = !bo;
            let v = if bo { -u } else { u };
            s = s * hh + v;
        }

        if bo {
            s = -s;
        }
        u = 0.5 * (u + 1.0);
        return Ok((u - s * h) * h.sqrt() * x + u * x.asinh());
    }

    w = 1.0 + f;
    if w == 0.0 {
        return Err("el3: 1 + px² cannot be zero.");
    }

    let p1 = if p == 0.0 { C::cb() / hh } else { p };
    s = s.abs();
    y = x_abs;
    g = if p == 1.0 { C::cb() } else { p1 - 1.0 };
    f = p1 - t;
    if f == 0.0 {
        f = C::cb() * t;
    }
    let am = 1.0 - t;
    let ap = 1.0 + e;
    r = p1 * h;
    fa = g / (f * p1);
    bo = fa > 0.0;
    fa = fa.abs();
    pz = (g * f).abs();
    de = pz.sqrt();
    q = p1.abs().sqrt();

    pm = pm.min(0.5);
    pm = p1 - pm;

    if pm >= 0.0 {
        u = (r * ap).sqrt();
        v = y * de;
        if g < 0.0 {
            v = -v;
        }
        d = 1.0 / q;
        c = 1.0;
    } else {
        u = (h * ap * pz).sqrt();
        ye = y * q;
        v = am * ye;
        q = -de / g;
        d = -am / de;
        c = 0.0;
        pz = ap - r;
    }

    if bo {
        r = v / u;
        z = 1.0;
        k = 1;
        if pm < 0.0 {
            h = y * (h / ap / fa).sqrt();
            h = 1.0 / h - h;
            z = h - r - r;
            r = 2.0 + r * h;
            if r == 0.0 {
                r = C::cb();
            }
            if z == 0.0 {
                z = h * C::cb();
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
        if p1 < 0.0 {
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
            ye = 0.0;
            u = (e + p1) / t;
            v = t / w;
            z = 1.0;
        }
        if s > 1.0 {
            h = u;
            u = v;
            v = h;
        }
    }

    y = 1.0 / y;
    e = s;
    nn = 1;
    t = 1.0;
    l = 0;
    mm = 0;

    for _ in 0..N_MAX_ITERATIONS {
        y = y - e / y;
        if y == 0.0 {
            y = e.sqrt() * C::cb();
        }
        f = c;
        c = d / q + c;
        g = e / q;
        d = f * g + d;
        d = d + d;
        q = g + q;
        g = t;
        t = s + t;
        nn = nn + nn;
        mm = mm + mm;
        if bo {
            if z < 0.0 {
                mm += k;
            }
            k = if r < 0.0 { -1 } else { 1 };
            h = e / (u * u + v * v);
            u = u * (1.0 + h);
            v = v * (1.0 - h);
        } else {
            r = u / v;
            h = z * r;
            z = h * z;
            hh = e / v;

            de = de / u;
            ye = ye * (h + 1.0 / h) + de * (1.0 + r);
            de = de * (u - hh);
            bk = ye.abs() < 1.0;

            // Removed crack function
        }
        if (g - s).abs() > C::ca() * g {
            if bo {
                g = (1.0 / r - r) * 0.5;
                hh = u + v * g;
                h = g * u - v;
                if hh == 0.0 {
                    hh = u * C::cb();
                }
                if h == 0.0 {
                    h = v * C::cb();
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
            if y < 0.0 {
                l += 1;
            }
            continue;
        }
        if y < 0.0 {
            l += 1;
        }
        e = (t / y).atan() + pi!() * T::from(l).unwrap();
        e = e * (c * t + d) / (t * (t + q));

        if bo {
            h = v / (t + u);
            z = 1.0 - r * h;
            h = r + h;
            if z == 0.0 {
                z = C::cb();
            }
            if z < 0.0 {
                mm += if h < 0.0 { -1 } else { 1 };
            }
            s = (h / z).atan() + T::from(mm).unwrap() * pi!();
        } else {
            s = if bk {
                ye.asinh()
            } else {
                z.ln() + T::from(mm).unwrap() * ln_2!()
            };
            s = s * 0.5;
        }
        e = (e + fa.sqrt() * s) / T::from(nn).unwrap();
        let ans = x.signum() * e;
        if ans.is_finite() {
            return Ok(ans);
        }
        check!(@nan, el3, [x, kc, p]);
        break;
    }

    let ans = ellippiinc(phi, n, m);
    if ans.is_ok() {
        return ans;
    }
    Err("el3: Failed to converge.")
}

const MAX_ND: usize = 50;

#[cfg(not(feature = "test_force_fail"))]
const N_MAX_ITERATIONS: usize = 10;

#[cfg(feature = "test_force_fail")]
const N_MAX_ITERATIONS: usize = 1;

#[cfg(not(feature = "test_force_fail"))]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{assert_close, compare_test_data_wolfram};

    #[test]
    fn test_el1() {
        compare_test_data_wolfram!("el1_data.csv", el1, 2, 5.0 * f64::EPSILON);
    }

    #[test]
    fn test_el1_references() {
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
    fn test_el1_special_cases() {
        use crate::bulirsch::cel1;
        use std::f64::{INFINITY, NAN};
        // x = 0: el1(0, kc) = 0
        assert_eq!(el1(0.0, 0.5).unwrap(), 0.0);
        // kc = 0: should return Err
        assert_eq!(el1(0.5, 0.0), Err("el1: kc cannot be zero."));
        // x = inf: el1(inf, kc) = cel1(kc)
        assert_eq!(el1(INFINITY, 0.5).unwrap(), cel1(0.5).unwrap());
        // kc = inf: el1(x, inf) = 0
        assert_eq!(el1(0.5, INFINITY).unwrap(), 0.0);
        // y = 0 branch in the loop
        assert_close!(el1(1.0, 1.0).unwrap(), 0.7853981633974483, 1e-15);
        // x = nan or kc = nan: should return Err
        assert_eq!(el1(NAN, 0.5), Err("el1: Arguments cannot be NAN."));
        assert_eq!(el1(0.5, NAN), Err("el1: Arguments cannot be NAN."));
    }

    #[test]
    fn test_el2() {
        compare_test_data_wolfram!("el2_data.csv", el2, 4, 100.0 * f64::EPSILON);
    }

    #[test]
    fn test_el2_special_cases() {
        use crate::bulirsch::cel2;
        use std::f64::{INFINITY, NAN};
        // x = 0: el2(0, kc, a, b) = 0
        assert_eq!(el2(0.0, 0.5, 1.0, 1.0).unwrap(), 0.0);
        // kc = 0: should return Err
        assert_eq!(el2(0.5, 0.0, 1.0, 1.0), Err("el2: kc cannot be zero."));
        // a = 0, b = 0: el2(x, kc, 0, 0) = 0
        assert_eq!(el2(0.5, 0.5, 0.0, 0.0).unwrap(), 0.0);
        // x = inf: el2(inf, kc, a, b) = cel2(kc, a, b)
        assert_eq!(
            el2(INFINITY, 0.5, 1.0, 1.0).unwrap(),
            cel2(0.5, 1.0, 1.0).unwrap()
        );
        // x = nan or kc = nan: should return Err
        assert_eq!(
            el2(NAN, 0.5, 1.0, 1.0),
            Err("el2: Arguments cannot be NAN.")
        );
        assert_eq!(
            el2(0.5, NAN, 1.0, 1.0),
            Err("el2: Arguments cannot be NAN.")
        );
    }

    #[test]
    fn test_el3() {
        compare_test_data_wolfram!("el3_data.csv", el3, 3, 3e-12);
        compare_test_data_wolfram!("el3_pv.csv", el3, 3, 5e-15);
    }

    #[test]
    fn test_el3_references() {
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
        test_reference(1.6, 1.50, 2.24999, 7.0057431688357934e-1);
        test_reference(1.6, 1e10, 1.20, 2.3734774669772208e-9);
        test_reference(-1.6, 1e10, 1.20, -2.3734774669772208e-9);
        test_reference(1.0, 0.31, 9.90e-2, 1.0903577921777398);

        // Use wolfram value at 16-digit precision instead since the values provided by
        // reference were calculated using 36-bit Matissa float (f32 and f64) making them inaccurate.
        test_reference(1.6, 1.90, -0.21, 1.473043981995436);
        test_reference(1.6, 1.90, -4.30, 0.2547695139719361);
        assert_close!(0.3950170978760504, el3(1.6, 1.01e1, -1.0e-5).unwrap(), 1e-9);
    }

    #[test]
    fn test_el3_special_cases() {
        use crate::cel;
        use std::f64::{INFINITY, NAN};
        // x = 0: el3(0, kc, p) = 0
        assert_eq!(el3(0.0, 0.5, 0.5).unwrap(), 0.0);
        // kc = 0: should return Err
        assert_eq!(el3(0.5, 0.0, 0.5), Err("el3: kc must not be zero."));

        let complete_el3 = el3(INFINITY, 0.5, 0.5).unwrap();
        // x = inf: el3(inf, kc, p) = cel(kc, p, 1, 1)
        assert_close!(complete_el3, cel(0.5, 0.5, 1.0, 1.0).unwrap(), 1e-15);
        // x = inf: el3(inf, kc, p) = Π(n, m)
        assert_close!(complete_el3, ellippi(0.5, 0.75).unwrap(), 1e-15);

        // kc = 1, p > 0: el3(x, 1, p) = atan(sqrt(p) * x) / sqrt(p)
        assert_eq!(
            el3(4.0, 1.0, 0.5).unwrap(),
            (0.5.sqrt() * 4.0).atan() / 0.5.sqrt()
        );
        // kc = 1, p <= 0: el3(x, 1, p) = (ln(1+vx) - ln(1-vx)) / (2v); v = sqrt(-p)
        assert_close!(el3(4.0, 1.0, -0.5).unwrap(), 5.0, 1e-15);
        // x = nan, kc = nan, or p = nan: should return Err
        assert_eq!(el3(NAN, 0.5, 0.5), Err("el3: Arguments cannot be NAN."));
        assert_eq!(el3(0.5, NAN, 0.5), Err("el3: Arguments cannot be NAN."));
        assert_eq!(el3(0.5, 0.5, NAN), Err("el3: Arguments cannot be NAN."));
    }
}

#[cfg(feature = "test_force_fail")]
crate::test_force_unreachable! {
    assert_eq!(el1(0.5, 0.5), Err("el1: Failed to converge."));
    assert_eq!(el2(0.5, 0.5, 0.5, 0.5), Err("el2: Failed to converge."));
    assert_eq!(el3(0.5, 0.5, 0.5), Err("el3: Failed to converge."));
}
