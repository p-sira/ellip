/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use num_traits::Float;

use crate::{
    bulirsch::{cel2::cel2_with_const, constants::BulirschConst, N_MAX_ITERATIONS},
    crate_util::{case, check, declare, let_mut},
    StrErr,
};

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
/// ![Bulirsch's Incomplete Elliptic Integral of the Second Kind](https://github.com/p-sira/ellip/blob/main/figures/el2.svg?raw=true)
///
/// [Interactive Plot](https://p-sira.github.io/ellippy/_static/figures/el2.html)
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
/// for [f64] and [f32]. The precision can be modified in the function [el2_with_const] (requires `unstable` feature flag).
///
/// # References
/// - Bulirsch, Roland. “Numerical Calculation of Elliptic Integrals and Elliptic Functions.” Numerische Mathematik 7, no. 1 (February 1, 1965): 78–90. <https://doi.org/10.1007/BF01397975>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
pub fn el2<T: Float>(x: T, kc: T, a: T, b: T) -> Result<T, StrErr> {
    if core::mem::size_of::<T>() <= 4 {
        el2_with_const::<T, f32>(x, kc, a, b)
    } else {
        el2_with_const::<T, f64>(x, kc, a, b)
    }
}

/// Computes [el2]. Control the precision using [BulirschConst].
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub fn el2_with_const<T: Float, C: BulirschConst<T>>(x: T, kc: T, a: T, b: T) -> Result<T, StrErr> {
    let ans = el2_unchecked::<T, C>(x, kc, a, b);
    if ans.is_finite() {
        return Ok(ans);
    }
    check!(@nan, el2, [x, kc, a, b]);
    check!(@zero, el2, [kc]);
    case!(x == T::zero(), T::zero());
    if x == inf!() {
        // phi = π/2
        return cel2_with_const::<T, C>(kc, a, b);
    }

    Err("el2: Failed to converge.")
}

/// Unsafe version of [el2].
/// <div class="warning">⚠️ Unstable feature. May subject to changes.</div>
///
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

#[cfg(not(feature = "test_force_fail"))]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::compare_test_data_wolfram;

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
}

#[cfg(feature = "test_force_fail")]
crate::test_force_unreachable! {
    use crate::bulirsch::constants::DefaultPrecision;
    assert_eq!(el2_with_const::<f64, DefaultPrecision>(0.5, 0.5, 0.5, 0.5), Err("el2: Failed to converge."));
}
