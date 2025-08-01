/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Boost Math, see LICENSE in this directory.
 */

// Original header from Boost Math
//  Copyright (c) 2006 Xiaogang Zhang
//  Copyright (c) 2006 John Maddock
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0.

use num_traits::Float;

use crate::crate_util::check;
use crate::legendre::ellippi::ellippi_vc;

use crate::{ellipeinc, ellipf, elliprc, elliprf, elliprj, StrErr};

/// Computes [incomplete elliptic integral of the third kind](https://dlmf.nist.gov/19.2.E7).
/// ```text
///                 φ                              
///                ⌠                 dϑ              
/// Π(φ, n, m)  =  ⎮ ──────────────────────────────────
///                ⎮   _____________                
///                ⌡ ╲╱ 1 - m sin²ϑ  ⋅ ( 1 - n sin²ϑ )
///               0              
/// ```
///
/// ## Parameters
/// - phi: amplitude angle (φ). φ ∈ ℝ.
/// - n: characteristic, n ∈ ℝ, n ≠ 1.
/// - m: elliptic parameter. m ∈ ℝ.
///
/// The elliptic modulus (k) is frequently used instead of the parameter (m), where k² = m.
/// The characteristic (n) is also sometimes expressed in term of α, where α² = n.
///
/// ## Domain
/// - Returns error when:
///   - m sin²φ > 1,
///   - n sin²φ = 1,
///   - or m ≥ 1 and φ is not a multiple of π/2.
/// - Returns the Cauchy principal value if n sin²φ > 1.
///
/// ## Graph
/// ![Incomplete Elliptic Integral of the Third Kind (3D Plot)](https://github.com/p-sira/ellip/blob/main/figures/ellippiinc_plot_3d.png?raw=true)
///
/// [Interactive 3D Plot](https://github.com/p-sira/ellip/blob/main/figures/ellippiinc_plot_3d.html)
///
/// ## Special Cases
/// - Π(0, n, m) = 0
/// - Π(φ, 0, 0) = φ
/// - Π(φ, 1, 0) = tan(φ)
/// - Π(φ, 0, m) = F(φ, m)
///
/// ## Notes
/// The [ellippiinc] (of the [crate]) is the circular or hyperbolic case of the elliptic integral
/// of the third kind, since n and m are real. It is called circular if `n (n - m) (n - 1)` is
/// negative and hyperbolic if it is positive.  
///
/// # Related Functions
/// - [ellippiinc](crate::ellippiinc)(φ, 0, m) = [ellipf](crate::ellipf)(φ, m)
/// - [ellippiinc](crate::ellippiinc)(φ, 1, m) = [sqrt(1 - m sin²(φ)) tan(φ) - [ellipe](crate::ellipe)(φ, m)] / (1 - m) + [ellipf](crate::ellipf)(φ, m)
/// - [ellippiinc](crate::ellippiinc)(φ, n, 1) = [sqrt(n) atanh(sqrt(n) sin(φ)) - log(sec(φ) + tan(φ))] / (n - 1)
///
/// With c = csc²φ,
/// - [ellippiinc](crate::ellippiinc)(φ, n, m) = n / 3 * [elliprj](crate::elliprj)(c - 1, c - m, c, c - n) + [ellipf](crate::ellipf)(φ, m)
///
/// With x = tan φ, p = 1 - n, and kc² = 1 - m,
/// - [ellippiinc](crate::ellippiinc)(φ, n, m) = [el3](crate::el3)(x, kc, p)
///
/// # Examples
/// ```
/// use ellip::{ellippiinc, util::assert_close};
/// use std::f64::consts::FRAC_PI_4;
///
/// assert_close(ellippiinc(FRAC_PI_4, 0.5, 0.5).unwrap(), 0.9190227391656969, 1e-15);
/// ```
///
/// # References
/// - Maddock, John, Paul Bristow, Hubert Holin, and Xiaogang Zhang. “Boost Math Library: Special Functions - Elliptic Integrals.” Accessed April 17, 2025. <https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/ellint.html>.
/// - Carlson, B. C. “DLMF: Chapter 19 Elliptic Integrals.” Accessed February 19, 2025. <https://dlmf.nist.gov/19>.
/// - Wolfram Research. “EllipticPi,” 2022. <https://reference.wolfram.com/language/ref/EllipticPi.html>.
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
pub fn ellippiinc<T: Float>(phi: T, n: T, m: T) -> Result<T, StrErr> {
    check!(@nan, ellippiinc, [phi, n, m]);
    ellippiinc_vc(phi, n, m, 1.0 - n)
}

#[inline]
#[numeric_literals::replace_float_literals(T::from(literal).unwrap())]
fn ellippiinc_vc<T: Float>(phi: T, n: T, m: T, nc: T) -> Result<T, StrErr> {
    // Note vc = 1-v presumably without cancellation error
    let sphi = phi.abs().sin();
    let sp2 = sphi * sphi;
    let mut result = 0.0;

    if m * sp2 > 1.0 {
        return Err("ellippiinc: m sin²φ must be smaller or equal to one.");
    }

    if n * sp2 == 1.0 {
        return Err("ellippiinc: n sin²φ must not equal one.");
    }

    // Special cases first:
    if n == 0.0 {
        // A&S 17.7.18 & 19
        return if m == 0.0 { Ok(phi) } else { ellipf(phi, m) };
    }

    if n == 1.0 {
        if m == 0.0 {
            return Ok(phi.tan());
        }

        // http://functions.wolfram.com/08.06.03.0008.01
        result = (1.0 - m * sp2).sqrt() * phi.tan() - ellipeinc(phi, m)?;
        result = result / (1.0 - m);
        result = result + ellipf(phi, m)?;
        return Ok(result);
    }

    if phi == pi_2!() {
        // Have to filter this case out before the next
        // special case, otherwise we might get an infinity from
        // tan(phi).
        // Also note that since we can't represent PI/2 exactly
        // this is a bit of a guess as to the users true intent...
        return ellippi_vc(n, m, nc);
    }

    // https://dlmf.nist.gov/19.6.E11
    if phi == 0.0 {
        return Ok(0.0);
    }

    if phi > pi_2!() || phi < 0.0 {
        // Carlson's algorithm works only for |phi| <= pi/2,
        // use the integrand's periodicity to normalize phi

        if phi.abs() > 1.0 / epsilon!() {
            // Invalid for v > 1, this case is caught above since v > 1 implies 1/v < sin^2(phi)
            debug_assert!(n <= 1.0);

            // Phi is so large that phi%pi is necessarily zero (or garbage),
            // just return the second part of the duplication formula:
            result = 2.0 * phi.abs() * ellippi_vc(n, m, nc)? / pi!();
        } else {
            let mut rphi = phi.abs() % pi_2!();
            let mut mm = ((phi.abs() - rphi) / pi_2!()).round();
            let mut sign = 1.0;

            if mm != 0.0 && m >= 1.0 {
                return Err("ellippiinc: The result is complex.");
            }

            if mm % 2.0 > 0.5 {
                mm = mm + 1.0;
                sign = -1.0;
                rphi = pi_2!() - rphi;
            }

            result = sign * ellippiinc_vc(rphi, n, m, nc)?;
            if mm > 0.0 && nc > 0.0 {
                result = result + mm * ellippi_vc(n, m, nc)?;
            }
        }
        return if phi < 0.0 { Ok(-result) } else { Ok(result) };
    }

    // https://reference.wolfram.com/language/ref/EllipticPi.html
    // Return the Cauchy principal value after normalizing phi
    if n * sp2 > 1.0 {
        let c = 1.0 / sp2; //csc2 phi
        let w2 = m / n;

        // This appears to have lower error in test dataset.
        // https://dlmf.nist.gov/19.7.E8
        return Ok((ellipf(phi, m)?
            + c.sqrt() * elliprc((c - 1.0) * (c - m), (c - n) * (c - w2))?)
            - ellippiinc(phi, w2, m)?);

        // https://dlmf.nist.gov/19.25.E16
        // return Ok(-third!() * w2 * elliprj(c - 1.0, c - m, c, c - w2)?
        //     + ((c - 1.0) * (c - m) / (n - 1.0) / (1.0 - w2)).sqrt()
        //         * elliprc(c * (n - 1.0) * (1.0 - w2), (n - c) * (c - w2))?);
    }

    if m == 0.0 {
        // A&S 17.7.20:
        if n < 1.0 {
            let vcr = nc.sqrt();
            return Ok((vcr * phi.tan()).atan() / vcr);
        } else {
            // v > 1:
            let vcr = (-nc).sqrt();
            let arg = vcr * phi.tan();
            return Ok((arg.ln_1p() - (-arg).ln_1p()) / (2.0 * vcr));
        }
    }

    if n < 0.0 && m <= 1.0 {
        //
        // If we don't shift to 0 <= v <= 1 we get
        // cancellation errors later on.  Use
        // A&S 17.7.15/16 to shift to v > 0.
        //
        // Mathematica simplifies the expressions
        // given in A&S as follows (with thanks to
        // Rocco Romeo for figuring these out!):
        //
        // V = (k2 - n)/(1 - n)
        // Assuming[(k2 >= 0 && k2 <= 1) && n < 0, FullSimplify[Sqrt[(1 - V)*(1 - k2 / V)] / Sqrt[((1 - n)*(1 - k2 / n))]]]
        // Result: ((-1 + k2) n) / ((-1 + n) (-k2 + n))
        //
        // Assuming[(k2 >= 0 && k2 <= 1) && n < 0, FullSimplify[k2 / (Sqrt[-n*(k2 - n) / (1 - n)] * Sqrt[(1 - n)*(1 - k2 / n)])]]
        // Result : k2 / (k2 - n)
        //
        // Assuming[(k2 >= 0 && k2 <= 1) && n < 0, FullSimplify[Sqrt[1 / ((1 - n)*(1 - k2 / n))]]]
        // Result : Sqrt[n / ((k2 - n) (-1 + n))]
        //
        let nn = (m - n) / (1.0 - n);
        let nm1 = (1.0 - m) / (1.0 - n);

        if nn > m {
            result = if nn.abs() < n.abs() {
                ellippiinc_vc(phi, nn, m, nm1)?
            } else {
                let c = 1.0 / sp2;
                nn / 3.0 * elliprj(c - 1.0, c - m, c, c - nn)? + ellipf(phi, m)?
            };
            result = result * n / (n - 1.0);
            result = result * (m - 1.0) / (n - m);
        }

        if m != 0.0 {
            let mut t = ellipf(phi, m)?;
            t = t * m / (m - n);
            result = result + t;
        }

        let mut p2 = -n * nn;
        if p2 <= min_val!() {
            p2 = (-n).sqrt() * nn.sqrt();
        } else {
            p2 = p2.sqrt();
        }

        let delta = (1.0 - m * sp2).sqrt();
        let t = n / ((m - n) * (n - 1.0));

        debug_assert!(t > min_val!());
        result = result + ((p2 / 2.0) * (2.0 * phi).sin() / delta).atan() * t.sqrt();

        // It should not be possible for t to be less than MIN value.
        // if t > min_val!() {
        //     result = result + ((p2 / 2.0) * (2.0 * phi).sin() / delta).atan() * t.sqrt();
        // } else {
        //     result = result
        //         + (((p2 / 2.0) * (2.0 * phi).sin() / delta).atan()
        //             * ((1.0 / (m - n)).abs()).sqrt()
        //             * ((n / (n - 1.0).abs()).sqrt()));
        // }

        return Ok(result);
    }

    if m == 1.0 {
        // See http://functions.wolfram.com/08.06.03.0013.01
        result = n.sqrt() * (n.sqrt() * phi.sin()).atanh() - (1.0 / phi.cos() + phi.tan()).ln();
        result = result / (n - 1.0);
        return Ok(result);
    }
    // disabled but retained for future reference: see below.
    //     if(v > 1)
    //    {
    //       //
    //       // If v > 1 we can use the identity in A&S 17.7.7/8
    //       // to shift to 0 <= v <= 1.  In contrast to previous
    //       // revisions of this header, this identity does now work
    //       // but appears not to produce better error rates in
    //       // practice.  Archived here for future reference...
    //       //
    //       T k2 = k * k;
    //       T N = k2 / v;
    //       T Nm1 = (v - k2) / v;
    //       T p1 = sqrt((-vc) * (1 - k2 / v));
    //       T delta = sqrt(1 - k2 * sphi * sphi);
    //       //
    //       // These next two terms have a large amount of cancellation
    //       // so it's not clear if this relation is useable even if
    //       // the issues with phi > pi/2 can be fixed:
    //       //
    //       result = -ellint_pi_imp(N, phi, k, Nm1, pol);
    //       result += ellint_f_imp(phi, k, pol);
    //       //
    //       // This log term gives the complex result when
    //       //     n > 1/sin^2(phi)
    //       // However that case is dealt with as an error above,
    //       // so we should always get a real result here:
    //       //
    //       result += log((delta + p1 * tan(phi)) / (delta - p1 * tan(phi))) / (2 * p1);
    //       return result;
    //    }

    //
    // Carlson's algorithm works only for |phi| <= pi/2,
    // by the time we get here phi should already have been
    // normalised above.
    //
    let cosp = phi.cos();
    let x = cosp * cosp;
    let t = sp2;
    let y = 1.0 - m * t;
    let z = 1.0;

    let p = if n * t < 0.5 { 1.0 - n * t } else { x + nc * t };

    let result = sphi * (elliprf(x, y, z)? + n * t * elliprj(x, y, z, p)? / 3.0);

    Ok(result)
}

#[cfg(not(feature = "reduce-iteration"))]
#[cfg(test)]
mod tests {
    use crate::compare_test_data_boost;

    use super::*;

    fn ellippiinc_k(inp: &[f64]) -> f64 {
        ellippiinc(inp[1], inp[0], inp[2] * inp[2]).unwrap()
    }

    #[test]
    fn test_ellippi() {
        compare_test_data_boost!("ellippiinc_data.txt", ellippiinc_k, 2.1e-15);
    }

    #[test]
    fn test_ellippi_large() {
        compare_test_data_boost!("ellippi3_large_data.txt", ellippiinc_k, 6e-15);
    }

    #[test]
    fn test_ellippiinc_err() {
        use std::f64::consts::FRAC_PI_2;
        // m sin²φ > 1
        assert!(ellippiinc(FRAC_PI_2, 0.5, 1.1).is_err());
        // n sin²φ = 1
        assert!(ellippiinc(FRAC_PI_2, 1.0, 0.5).is_err());
    }

    #[test]
    fn test_ellippiinc_special_cases() {
        use crate::ellippi;
        use std::f64::{
            consts::{FRAC_PI_2, PI},
            INFINITY, NAN, NEG_INFINITY,
        };
        // m * sin^2(phi) >= 1: should return Err
        assert!(ellippiinc(FRAC_PI_2, 0.5, 1.1).is_err());
        // n * sin^2(phi) = 1: should return Err
        assert!(ellippiinc(FRAC_PI_2, 1.0, 0.5).is_err());
        // Π(phi, 0, 0) = phi
        assert_eq!(ellippiinc(0.4, 0.0, 0.0).unwrap(), 0.4);
        // phi = 0: Π(0, n, m) = 0
        assert_eq!(ellippiinc(0.0, 0.5, 0.5).unwrap(), 0.0);
        // Π(phi, 1, 0) = tan(phi)
        assert_eq!(ellippiinc(0.2, 1.0, 0.0).unwrap(), 0.2.tan());
        // n = 0: Π(φ, 0, m) = F(φ, m)
        assert_eq!(
            ellippiinc(0.2, 0.0, 0.5).unwrap(),
            ellipf(0.2, 0.5).unwrap()
        );
        // n = 1: Π(φ, 1, m) = [sqrt(1-m sin2(phi)) tan(phi) - E(phi,m)]/(1-m) + F(phi,m)
        assert_eq!(
            ellippiinc(0.2, 1.0, 0.5).unwrap(),
            ((1.0 - 0.5 * 0.2.sin() * 0.2.sin()).sqrt() * 0.2.tan() - ellipeinc(0.2, 0.5).unwrap())
                / 0.5
                + ellipf(0.2, 0.5).unwrap()
        );
        // m = 1: Π(φ, n, 1) = [sqrt(n) atanh(sqrt(n) sin(φ)) - log(sec(φ)+tan(φ))]/(n-1)
        assert_eq!(
            ellippiinc(0.3, 1.5, 1.0).unwrap(),
            ((1.5.sqrt() * (1.5.sqrt() * 0.3.sin()).atanh()
                - (0.3.cos().recip() + 0.3.tan()).ln())
                / (1.5 - 1.0))
        );
        // phi > 1/epsilon: Π(φ, n, m) = 2 |phi| Π(n, m) / pi
        let large_phi_ans = 2.0 * 1e16 * ellippi(0.2, 0.5).unwrap() / PI;
        assert_eq!(ellippiinc(1e16, 0.2, 0.5).unwrap(), large_phi_ans);
        assert_eq!(ellippiinc(-1e16, 0.2, 0.5).unwrap(), -large_phi_ans);
        // phi = inf: Π(φ, n, m) = inf
        assert_eq!(ellippiinc(INFINITY, 0.2, 0.5).unwrap(), INFINITY);
        // phi = -inf: Π(φ, n, m) = -inf
        assert_eq!(ellippiinc(NEG_INFINITY, 0.2, 0.5).unwrap(), -INFINITY);
        // phi % pi/2 !=0, m >= 1: should return Err
        assert!(ellippiinc(4.14159, 0.5, 1.0).is_err());
        // phi = nan or n = nan or m = nan: should return Err
        assert!(ellippiinc(NAN, 0.5, 0.5).is_err());
        assert!(ellippiinc(0.5, NAN, 0.5).is_err());
        assert!(ellippiinc(0.5, 0.5, NAN).is_err());
    }
}
