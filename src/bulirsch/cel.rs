/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

#![allow(deprecated)]
use num_traits::Float;

use crate::{
    bulirsch::DefaultPrecision,
    crate_util::{case, check, declare},
    StrErr,
};

use super::_BulirschConst;

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
/// The precision of the function can be adjusted by overwriting the trait [super::BulirschConst].
/// The default is set according to the original literature by [Bulirsch](https://doi.org/10.1007/BF02165405) for [f64] and [f32].
///
/// ## Domain
/// - Returns error if kc = 0 or p = 0.
/// - Returns the Cauchy principal value for p < 0.
/// - Returns error if more than one arguments are infinite.
///
/// ## Graph
/// ![General Complete Elliptic Integral](https://github.com/p-sira/ellip/blob/main/figures/cel_plot.svg?raw=true)
///
/// [Interactive Plot](https://github.com/p-sira/ellip/blob/main/figures/cel_plot.html)
///
/// ## Special Cases
/// - cel(kc, p, 0, 0) = 0
/// - cel(kc, p, a, b) = 0 for |kc| = ∞
/// - cel(kc, p, a, b) = 0 for |p| = ∞
/// - cel(kc, p, 0, 0) = 0
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
/// # References
/// - Bulirsch, R. “Numerical Calculation of Elliptic Integrals and Elliptic Functions. III.” Numerische Mathematik 13, no. 4 (August 1, 1969): 305–15. <https://doi.org/10.1007/BF02165405>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
pub fn cel<T: Float>(kc: T, p: T, a: T, b: T) -> Result<T, StrErr> {
    _cel::<T, DefaultPrecision>(kc, p, a, b)
}

#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub fn _cel<T: Float, C: _BulirschConst<T>>(kc: T, p: T, a: T, b: T) -> Result<T, StrErr> {
    check!(@zero, cel, [kc, p]);

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
        aa = bb / pp + aa;
        g = e / pp;
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
/// The precision of the function can be adjusted by overwriting the trait [super::BulirschConst].
/// The default is set according to the original literature by [Bulirsch](https://doi.org/10.1007/BF02165405) for [f64] and [f32].
///
/// ## Domain
/// - Returns error if kc = 0.
///
/// ## Graph
/// ![Bulirsch's Complete Elliptic Integral of the First Kind](https://github.com/p-sira/ellip/blob/main/figures/cel1_plot.svg?raw=true)
///
/// [Interactive Plot](https://github.com/p-sira/ellip/blob/main/figures/cel1_plot.html)
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
/// # References
/// - Bulirsch, Roland. “Numerical Calculation of Elliptic Integrals and Elliptic Functions.” Numerische Mathematik 7, no. 1 (February 1, 1965): 78–90. <https://doi.org/10.1007/BF01397975>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
pub fn cel1<T: Float>(kc: T) -> Result<T, StrErr> {
    _cel1::<T, DefaultPrecision>(kc)
}

#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub fn _cel1<T: Float, C: _BulirschConst<T>>(kc: T) -> Result<T, StrErr> {
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
    if kc.is_infinite() {
        return Ok(0.0);
    }
    Err("cel1: Failed to converge.")
}

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
/// The precision of the function can be adjusted by overwriting the trait [super::BulirschConst].
/// The default is set according to the original literature by [Bulirsch](https://doi.org/10.1007/BF02165405) for [f64] and [f32].
///
/// ## Domain
/// - Returns error if kc = 0.
/// - Returns error if more than one arguments are infinite.
///
/// ## Graph
/// ![Bulirsch's Complete Elliptic Integral of the Second Kind](https://github.com/p-sira/ellip/blob/main/figures/cel2_plot.svg?raw=true)
///
/// [Interactive Plot](https://github.com/p-sira/ellip/blob/main/figures/cel2_plot.html)
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
/// # References
/// - Bulirsch, Roland. “Numerical Calculation of Elliptic Integrals and Elliptic Functions.” Numerische Mathematik 7, no. 1 (February 1, 1965): 78–90. <https://doi.org/10.1007/BF01397975>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
pub fn cel2<T: Float>(kc: T, a: T, b: T) -> Result<T, StrErr> {
    _cel2::<T, DefaultPrecision>(kc, a, b)
}

#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
#[inline]
pub fn _cel2<T: Float, C: _BulirschConst<T>>(kc: T, a: T, b: T) -> Result<T, StrErr> {
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

#[cfg(not(feature = "reduce-iteration"))]
const MAX_ITERATION: i16 = 10;
#[cfg(feature = "reduce-iteration")]
const MAX_ITERATION: i16 = 1;

#[cfg(not(feature = "reduce-iteration"))]
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
            assert_close!(ellipk, cel(kc, 1.0, 1.0, 1.0).unwrap(), 2e-12);
            assert_close!(ellipe, cel(kc, 1.0, 1.0, kc * kc).unwrap(), 1e-15);
            assert_close!(
                (ellipe - kc * kc * ellipk) / m,
                cel(kc, 1.0, 1.0, 0.0).unwrap(),
                8.5e-14
            );

            // Bulirsch, “Numerical Calculation of Elliptic Integrals and Elliptic Functions III”
            let n = 1.0 - p;
            let ellippi = ellippi(n, m).unwrap();

            assert_close!(
                (ellipk - ellipe) / m,
                cel(kc, 1.0, 0.0, 1.0).unwrap(),
                6e-12
            );
            // cel precision is very low for PI cases
            assert_close!(ellippi, cel(kc, p, 1.0, 1.0).unwrap(), 3.5e-12);
            assert_close!(
                (ellippi - ellipk) / (1.0 - p),
                cel(kc, p, 0.0, 1.0).unwrap(),
                3.5e-12
            );
        }

        let linsp_kc = [
            linspace(-1.0 + 1e-3, -1e-3, 100),
            linspace(1e-3, 1.0 - 1e-3, 100),
        ]
        .concat();
        let linsp_p = linspace(1e-3, 1.0 - 1e-3, 10);

        iproduct!(linsp_kc, linsp_p).for_each(|(kc, p)| _test(kc, p));

        // Data from Bulirsch, “Numerical Calculation of Elliptic Integrals and Elliptic Functions III”
        assert_close!(cel(1e-1, 4.1, 1.2, 1.1).unwrap(), 1.5464442694017956, 5e-16);
        assert_close!(
            cel(1e-1, -4.1, 1.2, 1.1).unwrap(),
            -6.7687378198360556e-1,
            5e-16
        );
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

#[cfg(feature = "reduce-iteration")]
crate::test_force_unreachable! {
            assert_eq!(cel(1e300, 0.2, 0.5, 0.5), Err("cel: Arguments cannot be NAN."));
        assert_eq!(cel1(1e300), Err("cel1: Arguments cannot be NAN."));
        assert_eq!(cel2(1e300, 0.5, 0.5), Err("cel2: Arguments cannot be NAN."));
}
