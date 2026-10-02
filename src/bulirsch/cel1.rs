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

/// Computes [complete elliptic integral of the first kind in Bulirsch's form](https://link.springer.com/article/10.1007/bf01397975).
/// ```text
///               π/2                                                   
///              ⌠               dϑ              
/// cel1(kc)  =  ⎮  ────────────────────────────
///              ⎮      ______________________
///              ⌡   ╲╱ cos²(ϑ) + kc² sin²(ϑ)    
///             0                                                   
/// ```
///
/// ## Parameters
/// - kc: complementary modulus. kc ∈ ℝ, kc ≠ 0.
///
/// ## Domain
/// - Returns error if kc = 0.
///
/// ## Graph
/// ![Bulirsch's Complete Elliptic Integral of the First Kind](https://github.com/p-sira/ellip/blob/main/figures/cel1.svg?raw=true)
///
/// [Interactive Plot](https://p-sira.github.io/ellippy/_static/figures/cel1.html)
///
/// ## Special Cases
/// - cel1(kc) = 0 for |kc| = ∞
///
/// # Related Functions
/// With kc² = 1 - m,
/// - [ellipk](crate::ellipk)(m) = [cel](crate::cel)(kc, 1, 1, 1) = [cel1](crate::cel1)(kc)
///
/// # Examples
/// ```
/// use ellip::{cel1, util::assert_close};
///
/// assert_close(cel1(0.5).unwrap(), 2.1565156474996434, 1e-15);
/// ```
///
///  # Notes
/// The default precision of the function is set according to the original literature by [Bulirsch](https://doi.org/10.1007/BF02165405)
/// for [f64] and [f32]. The precision can be modified in the function [cel1_with_const] (requires `unstable` feature flag).
///
/// # References
/// - Bulirsch, Roland. “Numerical Calculation of Elliptic Integrals and Elliptic Functions.” Numerische Mathematik 7, no. 1 (February 1, 1965): 78–90. <https://doi.org/10.1007/BF01397975>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
pub fn cel1<T: Float>(kc: T) -> Result<T, StrErr> {
    if core::mem::size_of::<T>() <= 4 {
        cel1_with_const::<T, f32>(kc)
    } else {
        cel1_with_const::<T, f64>(kc)
    }
}

/// Computes [cel1]. Control the precision using [BulirschConst].
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub fn cel1_with_const<T: Float, C: BulirschConst<T>>(kc: T) -> Result<T, StrErr> {
    declare!(mut [kc = kc.abs(), m = T::one(), ans = T::nan(), h]);
    for _ in 0..MAX_ITERATION {
        h = m;
        m = kc + m;

        if (h - kc).abs() > C::ca() * h {
            kc = (h * kc).sqrt();
            m = m / 2.0;
            continue;
        }

        ans = pi!() / m;
        break;
    }

    if ans.is_finite() {
        return Ok(ans);
    }
    check!(@nan, cel1, [kc]);
    check!(@zero, cel1, [kc]);
    Err("cel1: Failed to converge.")
}

#[cfg(not(feature = "test_force_fail"))]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{assert_close, ellipk, test_util::linspace};

    #[test]
    fn test_cel1() {
        fn test_kc(kc: f64) {
            let m = 1.0 - kc * kc;
            assert_close!(ellipk(m).unwrap(), cel1(kc).unwrap(), 2e-12);
        }

        let linsp_neg = linspace(-1.0, -1e-3, 100);
        linsp_neg.iter().for_each(|kc| test_kc(*kc));
        let linsp_pos = linspace(1e-3, 1.0, 100);
        linsp_pos.iter().for_each(|kc| test_kc(*kc));
    }

    #[test]
    fn test_cel1_special_cases() {
        use std::f64::{INFINITY, NAN, NEG_INFINITY};
        // kc = 0: should return Err
        assert_eq!(cel1(0.0), Err("cel1: kc cannot be zero."));
        // kc = inf: cel1(inf) = 0
        assert_eq!(cel1(INFINITY).unwrap(), 0.0);
        // kc = -inf: cel1(-inf) = 0
        assert_eq!(cel1(NEG_INFINITY).unwrap(), 0.0);
        // kc = NaN: should return Err
        assert_eq!(cel1(NAN), Err("cel1: Arguments cannot be NAN."));
    }
}

#[cfg(feature = "test_force_fail")]
crate::test_force_unreachable! {
    assert_eq!(cel1(1e300), Err("cel1: Failed to converge."));
}
