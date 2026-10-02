/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use num_traits::Float;

use crate::{
    bulirsch::{constants::BulirschConst, MAX_ITERATION},
    crate_util::{check, declare},
    StrErr,
};

/// Computes [complete elliptic integral of the second kind in Bulirsch's form](https://link.springer.com/article/10.1007/bf01397975).
/// ```text
///                     π/2                           
///                    ⌠              a + b tan²(ϑ)
/// cel2(kc, a, b)  =  ⎮  ──────────────────────────────────── ⋅ dϑ
///                    ⎮     ________________________________
///                    ⌡  ╲╱ (1 + tan²(ϑ)) (1 + kc² tan²(ϑ))    
///                    0                                                   
/// where kc ≠ 0
/// ```
///
/// ## Parameters
/// - kc: complementary modulus. kc ∈ ℝ, kc ≠ 0.
/// - a ∈ ℝ
/// - b ∈ ℝ
///
/// ## Domain
/// - Returns error if kc = 0.
/// - Returns error if more than one arguments are infinite.
///
/// ## Graph
/// ![Bulirsch's Complete Elliptic Integral of the Second Kind](https://github.com/p-sira/ellip/blob/main/figures/cel2.svg?raw=true)
///
/// [Interactive Plot](https://p-sira.github.io/ellippy/_static/figures/cel2.html)
///
/// ## Special Cases
/// - cel2(kc, 0, 0) = 0
/// - cel(kc, a, b) = 0 for |kc| = ∞
/// - cel(kc, a, b) = sign(a) ∞ for |a| = ∞
/// - cel(kc, a, b) = sign(b) ∞ for |b| = ∞
///
/// # Related Functions
/// - [cel2](crate::cel2)(kc, a, b) = [cel](crate::cel)(kc, 1, a, b)
///
/// With kc² = 1 - m,
/// - [ellipe](crate::ellipe)(m) = [cel](crate::cel)(kc, 1, 1, kc²) = [cel2](crate::cel2)(kc, 1, kc²)
///
/// # Examples
/// ```
/// use ellip::{cel2, util::assert_close};
///
/// assert_close(cel2(0.5, 1.0, 1.0).unwrap(), 2.1565156474996434, 1e-15);
/// ```
///
/// # Notes
/// The default precision of the function is set according to the original literature by [Bulirsch](https://doi.org/10.1007/BF02165405)
/// for [f64] and [f32]. The precision can be modified in the function [cel2_with_const] (requires `unstable` feature flag).
///
/// # References
/// - Bulirsch, Roland. “Numerical Calculation of Elliptic Integrals and Elliptic Functions.” Numerische Mathematik 7, no. 1 (February 1, 1965): 78–90. <https://doi.org/10.1007/BF01397975>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
pub fn cel2<T: Float>(kc: T, a: T, b: T) -> Result<T, StrErr> {
    if core::mem::size_of::<T>() <= 4 {
        cel2_with_const::<T, f32>(kc, a, b)
    } else {
        cel2_with_const::<T, f64>(kc, a, b)
    }
}

/// Computes [cel2]. Control the precision using [BulirschConst].
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub fn cel2_with_const<T: Float, C: BulirschConst<T>>(kc: T, a: T, b: T) -> Result<T, StrErr> {
    declare!(mut [kc = kc.abs(), aa = a, bb = b, m = T::one(), c = aa, ans = T::nan(), m0]);
    aa = bb + aa;

    for _ in 0..MAX_ITERATION {
        bb = (c * kc + bb) * 2.0;
        c = aa;
        m0 = m;
        m = kc + m;
        aa = bb / m + aa;

        if (m0 - kc).abs() > C::ca() * m0 {
            kc = (kc * m0).sqrt() * 2.0;
            continue;
        }

        ans = pi!() / 4.0 * aa / m;
        break;
    }

    if ans.is_finite() {
        return Ok(ans);
    }
    check!(@nan, cel2, [kc, a, b]);
    check!(@zero, cel2, [kc]);
    check!(@multi, cel2, "infinite", is_infinite, [kc, a, b]);
    if kc.is_infinite() {
        return Ok(0.0);
    }
    if a.is_infinite() {
        return Ok(a.signum() * inf!());
    }
    if b.is_infinite() {
        return Ok(b.signum() * inf!());
    }
    Err("cel2: Failed to converge.")
}

#[cfg(not(feature = "test_force_fail"))]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{assert_close, ellipe, test_util::linspace};

    #[test]
    fn test_cel2() {
        fn test_kc(kc: f64) {
            let m = 1.0 - kc * kc;
            assert_close!(ellipe(m).unwrap(), cel2(kc, 1.0, kc * kc).unwrap(), 7.7e-16);
        }

        let linsp_neg = linspace(-1.0, -1e-3, 100);
        linsp_neg.iter().for_each(|kc| test_kc(*kc));
        let linsp_pos = linspace(1e-3, 1.0, 100);
        linsp_pos.iter().for_each(|kc| test_kc(*kc));
    }

    #[test]
    fn test_cel2_special_cases() {
        use std::f64::{INFINITY, NAN, NEG_INFINITY};
        // kc = 0: should return Err
        assert_eq!(cel2(0.0, 1.0, 1.0), Err("cel2: kc cannot be zero."));
        // kc = inf: cel2(inf, 1, 1) = 0
        assert_eq!(cel2(INFINITY, 1.0, 1.0).unwrap(), 0.0);
        // kc = -inf: cel2(-inf, 1, 1) = 0
        assert_eq!(cel2(NEG_INFINITY, 1.0, 1.0).unwrap(), 0.0);
        // a = 0 and b = 0: cel2(kc, a, b) = 0
        assert_eq!(cel2(0.5, 0.0, 0.0).unwrap(), 0.0);
        // a = inf: cel2(kc, inf, b) = inf
        assert_eq!(cel2(0.5, INFINITY, 1.0).unwrap(), INFINITY);
        // b = inf: cel2(kc, a, inf) = inf
        assert_eq!(cel2(0.5, 1.0, INFINITY).unwrap(), INFINITY);
        // a = inf: cel2(kc, -inf, b) = -inf
        assert_eq!(cel2(0.5, NEG_INFINITY, 1.0).unwrap(), NEG_INFINITY);
        // b = inf: cel2(kc, a, -inf) = -inf
        assert_eq!(cel2(0.5, 1.0, NEG_INFINITY).unwrap(), NEG_INFINITY);
        // NANs: should return Err
        assert_eq!(cel2(NAN, 1.0, 1.0), Err("cel2: Arguments cannot be NAN."));
        assert_eq!(cel2(0.5, NAN, 1.0), Err("cel2: Arguments cannot be NAN."));
        assert_eq!(cel2(0.5, 1.0, NAN), Err("cel2: Arguments cannot be NAN."));
    }
}

#[cfg(feature = "test_force_fail")]
crate::test_force_unreachable! {
    assert_eq!(cel2(1e300, 0.5, 0.5), Err("cel2: Failed to converge."));
}
