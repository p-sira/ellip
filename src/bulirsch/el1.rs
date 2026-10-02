/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use num_traits::Float;

use crate::{
    bulirsch::{constants::BulirschConst, N_MAX_ITERATIONS},
    crate_util::{case, check, declare},
    StrErr,
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
/// ![Bulirsch's Incomplete Elliptic Integral of the First Kind](https://github.com/p-sira/ellip/blob/main/figures/el1.svg?raw=true)
///
/// [Interactive Plot](https://p-sira.github.io/ellippy/_static/figures/el1.html)
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
/// for [f64] and [f32]. The precision can be modified in the function [el1_with_const] (requires `unstable` feature flag).
///
/// # References
/// - Bulirsch, Roland. “Numerical Calculation of Elliptic Integrals and Elliptic Functions.” Numerische Mathematik 7, no. 1 (February 1, 1965): 78–90. <https://doi.org/10.1007/BF01397975>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
pub fn el1<T: Float>(x: T, kc: T) -> Result<T, StrErr> {
    if core::mem::size_of::<T>() <= 4 {
        el1_with_const::<T, f32>(x, kc)
    } else {
        el1_with_const::<T, f64>(x, kc)
    }
}

/// Computes [el1]. Control the precision using [BulirschConst].
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub fn el1_with_const<T: Float, C: BulirschConst<T>>(x: T, kc: T) -> Result<T, StrErr> {
    if x == 0.0 {
        check!(@nan, el1, [x, kc]);
        check!(@zero, el1, [kc]);
        return Ok(0.0);
    }
    let ans = el1_unchecked::<T, C>(x, kc);
    if ans.is_finite() {
        return Ok(ans);
    }
    check!(@nan, el1, [x, kc]);
    check!(@zero, el1, [kc]);
    case!(kc == inf!(), T::zero());
    Err("el1: Failed to converge.")
}

/// Unsafe version of [el1].
/// <div class="warning">⚠️ Unstable feature. May subject to changes.</div>
///
/// Undefined behavior with invalid arguments and edge cases.
/// # Known Invalid Cases
/// - kc = 0
/// - kc = ∞
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub fn el1_unchecked<T: Float, C: BulirschConst<T>>(x: T, kc: T) -> T {
    if x == 0.0 {
        return 0.0;
    }
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
        assert_close!(el1(1.0, 1.0).unwrap(), std::f64::consts::FRAC_PI_4, 1e-15);
        // x = nan or kc = nan: should return Err
        assert_eq!(el1(NAN, 0.5), Err("el1: Arguments cannot be NAN."));
        assert_eq!(el1(0.5, NAN), Err("el1: Arguments cannot be NAN."));
    }
}

#[cfg(feature = "test_force_fail")]
crate::test_force_unreachable! {
    use crate::bulirsch::constants::DefaultPrecision;
    assert_eq!(el1_with_const::<f64, DefaultPrecision>(0.5, 0.5), Err("el1: Failed to converge."));
}
