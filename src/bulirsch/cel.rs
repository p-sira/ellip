/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use num_traits::Float;

use crate::{
    bulirsch::{
        cel1::cel1_with_const, cel2::cel2_with_const, constants::BulirschConst, MAX_ITERATION,
    },
    crate_util::{case, check, declare},
    StrErr,
};

/// Computes [complete elliptic integral in Bulirsch form](https://dlmf.nist.gov/19.2#iii).
/// ```text
///                       π/2                                                   
///                      ⌠            a cos²(ϑ) + b sin²(ϑ)               
/// cel(kc, p, a, b)  =  ⎮ ─────────────────────────────────────────────── ⋅ dϑ
///                      ⎮                         ______________________
///                      ⌡ (cos²(ϑ) + p sin²(ϑ)) ╲╱ cos²(ϑ) + kc² sin²(ϑ)    
///                     0                                                   
/// ```
///
/// ## Parameters
/// - kc: complementary modulus. kc ∈ ℝ, kc ≠ 0.
/// - p ∈ ℝ, p ≠ 0
/// - a ∈ ℝ
/// - b ∈ ℝ
///
/// ## Domain
/// - Returns error if kc = 0 or p = 0.
/// - Returns the Cauchy principal value for p < 0.
/// - Returns error if more than one arguments are infinite.
///
/// ## Graph
/// ![General Complete Elliptic Integral](https://github.com/p-sira/ellip/blob/main/figures/cel.svg?raw=true)
///
/// [Interactive Plot](https://p-sira.github.io/ellippy/_static/figures/cel.html)
///
/// ## Special Cases
/// - cel(kc, p, 0, 0) = 0
/// - cel(kc, p, a, b) = 0 for |kc| = ∞
/// - cel(kc, p, a, b) = 0 for |p| = ∞
/// - cel(kc, p, a, b) = sign(a) ∞ for |a| = ∞
/// - cel(kc, p, a, b) = sign(b) ∞ for |b| = ∞
///
/// # Related Functions
/// With kc² = 1 - m and p = 1 - n,
/// - [ellipk](crate::ellipk)(m) = [cel](crate::cel)(kc, 1, 1, 1) = [cel1](crate::cel1)(kc)
/// - [ellipe](crate::ellipe)(m) = [cel](crate::cel)(kc, 1, 1, kc²) = [cel2](crate::cel2)(kc, 1, kc²)
/// - [ellipd](crate::ellipd)(m) = [cel](crate::cel)(kc, 1, 0, 1)
/// - [ellippi](crate::ellippi)(n, m) = [cel](crate::cel)(kc, p, 1, 1)
///
/// # Examples
/// ```
/// use ellip::{cel, util::assert_close};
///
/// assert_close(cel(0.5, 1.0, 1.0, 1.0).unwrap(), 2.1565156474996434, 1e-15);
/// ```
///
/// # Notes
/// The default precision of the function is set according to the original literature by [Bulirsch](https://doi.org/10.1007/BF02165405)
/// for [f64] and [f32]. The precision can be modified in the function [cel_with_const] (requires `unstable` feature flag).
///
/// # References
/// - Bulirsch, R. “Numerical Calculation of Elliptic Integrals and Elliptic Functions. III.” Numerische Mathematik 13, no. 4 (August 1, 1969): 305–15. <https://doi.org/10.1007/BF02165405>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
pub fn cel<T: Float>(kc: T, p: T, a: T, b: T) -> Result<T, StrErr> {
    if core::mem::size_of::<T>() <= 4 {
        cel_with_const::<T, f32>(kc, p, a, b)
    } else {
        cel_with_const::<T, f64>(kc, p, a, b)
    }
}

/// Computes [cel]. Control the precision using [BulirschConst].
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub fn cel_with_const<T: Float, C: BulirschConst<T>>(kc: T, p: T, a: T, b: T) -> Result<T, StrErr> {
    check!(@zero, cel, [kc, p]);

    if p == 1.0 {
        let res = if a == 1.0 && b == 1.0 {
            cel1_with_const::<T, C>(kc)
        } else {
            cel2_with_const::<T, C>(kc, a, b)
        };
        if let Ok(v) = res {
            return Ok(v);
        }
    }

    let mut kc = kc.abs();
    declare!(mut [pp = p, aa = a, bb = b, f, q, g]);

    let mut e = kc;
    let mut m = 1.0;

    if pp > 0.0 {
        pp = pp.sqrt();
        bb = bb / pp;
    } else {
        f = kc * kc;
        q = 1.0 - f;
        g = 1.0 - pp;
        f = f - pp;
        q = (bb - aa * pp) * q;
        pp = (f / g).sqrt();
        aa = (aa - bb) / g;
        bb = -q / (g * g * pp) + aa * pp;
    }

    let mut ans = T::nan();
    for _ in 0..MAX_ITERATION {
        f = aa;
        let inv_pp = 1.0 / pp;
        aa = bb * inv_pp + aa;
        g = e * inv_pp;
        bb = 2.0 * (f * g + bb);
        pp = g + pp;
        g = m;
        m = kc + m;

        if (g - kc).abs() > g * C::ca() {
            kc = 2.0 * e.sqrt();
            e = kc * m;
            continue;
        }

        ans = pi_2!() * (aa * m + bb) / (m * (m + pp));
        break;
    }

    if ans.is_finite() {
        return Ok(ans);
    }
    check!(@nan, cel, [kc, p, a, b]);
    check!(@multi, cel, "infinite", is_infinite, [kc, p, a, b]);
    case!(@any [kc.abs(), p.abs()] == inf!(), T::zero());
    if a.is_infinite() {
        return Ok(a.signum() * inf!());
    }
    if b.is_infinite() {
        return Ok(b.signum() * inf!());
    }
    Err("cel: Failed to converge.")
}

#[cfg(not(feature = "test_force_fail"))]
#[cfg(test)]
mod tests {
    use itertools::iproduct;

    use super::*;
    use crate::{assert_close, ellipe, ellipk, ellippi, test_util::linspace};

    /// Test using relationship with Legendre form.
    /// Reference: https://dlmf.nist.gov/19.2#iii
    #[test]
    fn test_cel() {
        fn _test(kc: f64, p: f64) {
            let m = 1.0 - kc * kc;
            let ellipk = ellipk(m).unwrap();
            let ellipe = ellipe(m).unwrap();

            // cel precision is low for K cases
            assert_close! {ellipk, cel(kc, 1.0, 1.0, 1.0).unwrap(), 2e-12};
            assert_close! {ellipe, cel(kc, 1.0, 1.0, kc * kc).unwrap(), 1e-15};
            assert_close! {(ellipe - kc * kc * ellipk) / m, cel(kc, 1.0, 1.0, 0.0).unwrap(), 8.5e-14};

            // Bulirsch, “Numerical Calculation of Elliptic Integrals and Elliptic Functions III”
            let n = 1.0 - p;
            let ellippi = ellippi(n, m).unwrap();

            assert_close! {(ellipk - ellipe) / m, cel(kc, 1.0, 0.0, 1.0).unwrap(), 6e-12};
            // cel precision is very low for PI cases
            assert_close! {ellippi, cel(kc, p, 1.0, 1.0).unwrap(), 3.5e-12};
            assert_close! {(ellippi - ellipk) / (1.0 - p), cel(kc, p, 0.0, 1.0).unwrap(), 3.5e-12};
        }

        let linsp_kc = [
            linspace(-1.0 + 1e-3, -1e-3, 100),
            linspace(1e-3, 1.0 - 1e-3, 100),
        ]
        .concat();
        let linsp_p = linspace(1e-3, 1.0 - 1e-3, 10);

        iproduct!(linsp_kc, linsp_p).for_each(|(kc, p)| _test(kc, p));

        // Data from Bulirsch, “Numerical Calculation of Elliptic Integrals and Elliptic Functions III”
        assert_close! {cel(1e-1, 4.1, 1.2, 1.1).unwrap(), 1.5464442694017956, 5e-16};
        assert_close! {cel(1e-1, -4.1, 1.2, 1.1).unwrap(), -6.7687378198360556e-1, 5e-16};
    }

    #[test]
    fn test_cel_special_cases() {
        use std::f64::{INFINITY, NAN, NEG_INFINITY};
        // kc = 0: should return Err
        assert_eq!(cel(0.0, 1.0, 1.0, 1.0), Err("cel: kc cannot be zero."));
        // p = 0: should return Err
        assert_eq!(cel(0.5, 0.0, 1.0, 1.0), Err("cel: p cannot be zero."));
        // a = 0, b = 0: cel(kc, p, 0, 0) = 0
        assert_eq!(cel(0.5, 1.0, 0.0, 0.0).unwrap(), 0.0);
        // kc = inf: cel(inf, p, a, b) = 0
        assert_eq!(cel(INFINITY, 0.5, 1.0, 1.0).unwrap(), 0.0);
        // kc = -inf: cel(-inf, kc, 1.0, 1.0) = 0
        assert_eq!(cel(NEG_INFINITY, 0.5, 1.0, 1.0).unwrap(), 0.0);
        // p = inf: cel(kc, inf, a, b) = 0
        assert_eq!(cel(0.5, INFINITY, 1.0, 1.0).unwrap(), 0.0);
        // p = -inf: cel(kc, -inf, a, b) = 0
        assert_eq!(cel(0.5, NEG_INFINITY, 1.0, 1.0).unwrap(), 0.0);
        // a = inf: cel(kc, p, inf, b) = inf
        assert_eq!(cel(0.5, 1.0, INFINITY, 1.0).unwrap(), INFINITY);
        // b = inf: cel(kc, p, a, inf) = inf
        assert_eq!(cel(0.5, 1.0, 1.0, INFINITY).unwrap(), INFINITY);
        // a = inf: cel(kc, p, -inf, b) = -inf
        assert_eq!(cel(0.5, 1.0, NEG_INFINITY, 1.0).unwrap(), NEG_INFINITY);
        // b = inf: cel(kc, p, a, -inf) = -inf
        assert_eq!(cel(0.5, 1.0, 1.0, NEG_INFINITY).unwrap(), NEG_INFINITY);
        // NANs: should return Err
        assert_eq!(
            cel(NAN, 1.0, 1.0, 1.0),
            Err("cel: Arguments cannot be NAN.")
        );
        assert_eq!(
            cel(0.5, NAN, 1.0, 1.0),
            Err("cel: Arguments cannot be NAN.")
        );
        assert_eq!(
            cel(0.5, 1.0, NAN, 1.0),
            Err("cel: Arguments cannot be NAN.")
        );
        assert_eq!(
            cel(0.5, 1.0, 1.0, NAN),
            Err("cel: Arguments cannot be NAN.")
        );
    }
}

#[cfg(feature = "test_force_fail")]
crate::test_force_unreachable! {
    assert_eq!(cel(1e300, 0.2, 0.5, 0.5), Err("cel: Failed to converge."));
}
