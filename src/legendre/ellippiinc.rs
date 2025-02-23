/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Boost Math, see LICENSE in this directory.
 */

// Original header from Boost Math
//  Copyright (c) 2006 Xiaogang Zhang
//  Copyright (c) 2006 John Maddock
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  History:
//  XZ wrote the original of this file as part of the Google
//  Summer of Code 2006.  JM modified it to fit into the
//  Boost.Math conceptual framework better, and to correctly
//  handle the various corner cases.
//

use num_traits::Float;

use crate::legendre::ellippi::ellippi_vc;

use crate::{ellipeinc, ellipf, elliprf, elliprj};

/// Compute [incomplete elliptic integral of the third kind](https://dlmf.nist.gov/19.2.E7).
/// ```text
///                 φ                              
///                ⌠                 dϑ              
/// Π(φ, n, m)  =  ⎮ ──────────────────────────────────
///                ⎮   _____________                
///                ⌡ ╲╱ 1 - m sin²ϑ  ⋅ ( 1 - n sin²ϑ )
///               0              
/// where 0 ≤ m sin²φ ≤ 1               
/// ```
///
/// Note that some mathematical references use the parameter k and α for the function,
/// where k² = m, α² = n.
pub fn ellippiinc<T: Float>(phi: T, n: T, m: T) -> Result<T, &'static str> {
    ellippiinc_vc(phi, n, m, one!() - n)
}

#[inline]
fn ellippiinc_vc<T: Float>(phi: T, n: T, m: T, vc: T) -> Result<T, &'static str> {
    // Note vc = 1-v presumably without cancellation error
    let sphi = phi.abs().sin();
    let sp2 = sphi * sphi;
    let mut result = zero!();

    if m * sp2 > one!() || m < zero!() {
        return Err("ellippiinc: The argument must satisfy 0 ≤ m sin²φ ≤ 1.");
    }

    // Special cases first:
    if n == zero!() {
        // A&S 17.7.18 & 19
        return if m == zero!() {
            Ok(phi)
        } else {
            ellipf(phi, m)
        };
    }

    if n > zero!() && (one!() / n) < sp2 {
        // Complex result is a domain error:
        return Err("ellippiinc: result is complex for v > 1 / sin^2(phi)");
    }

    if n == one!() {
        if m == zero!() {
            return Ok(phi.tan());
        }

        // http://functions.wolfram.com/08.06.03.0008.01
        result = (one!() - m * sp2).sqrt() * phi.tan() - ellipeinc(phi, m)?;
        result = result / (one!() - m);
        result = result + ellipf(phi, m)?;
        return Ok(result);
    }

    if phi == pi_2!() {
        // Have to filter this case out before the next
        // special case, otherwise we might get an infinity from
        // tan(phi).
        // Also note that since we can't represent PI/2 exactly
        // this is a bit of a guess as to the users true intent...
        return ellippi_vc(n, m, vc);
    }

    if phi > pi_2!() || phi < zero!() {
        // Carlson's algorithm works only for |phi| <= pi/2,
        // use the integrand's periodicity to normalize phi

        if phi.abs() > one!() / epsilon!() {
            // Invalid for v > 1, this case is caught above since v > 1 implies 1/v < sin^2(phi)
            debug_assert!(n <= one!());

            // Phi is so large that phi%pi is necessarily zero (or garbage),
            // just return the second part of the duplication formula:
            result = two!() * phi.abs() * ellippi_vc(n, m, vc)? / pi!();
        } else {
            let mut rphi = phi.abs() % pi_2!();
            let mut mm = ((phi.abs() - rphi) / pi_2!()).round();
            let mut sign = one!();

            if mm != zero!() && m >= one!() {
                return Err("ellippiinc: The result is complex.");
            }

            if mm % two!() > half!() {
                mm = mm + one!();
                sign = -one!();
                rphi = pi_2!() - rphi;
            }

            result = sign * ellippiinc_vc(rphi, n, m, vc)?;
            if mm > zero!() && vc > zero!() {
                result = result + mm * ellippi_vc(n, m, vc)?;
            }
        }
        return if phi < zero!() {
            Ok(-result)
        } else {
            Ok(result)
        };
    }

    if m == zero!() {
        // A&S 17.7.20:
        if n < one!() {
            let vcr = vc.sqrt();
            return Ok((vcr * phi.tan()).atan() / vcr);
        } else {
            // v > 1:
            let vcr = (-vc).sqrt();
            let arg = vcr * phi.tan();
            return Ok((arg.ln_1p() - (-arg).ln_1p()) / (two!() * vcr));
        }
    }

    if n < zero!() && m <= one!() {
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
        let nn = (m - n) / (one!() - n);
        let nm1 = (one!() - m) / (one!() - n);
        let mut p2 = -n * nn;

        if p2 <= min_val!() {
            p2 = (-n).sqrt() * nn.sqrt();
        } else {
            p2 = p2.sqrt();
        }

        let delta = (one!() - m * sp2).sqrt();
        if nn > m {
            result = ellippiinc_vc(phi, nn, m, nm1)?;
            result = result * n / (n - one!());
            result = result * (m - one!()) / (n - m);
        }

        if m != zero!() {
            let mut t = ellipf(phi, m)?;
            t = t * m / (m - n);
            result = result + t;
        }
        let t = n / ((m - n) * (n - one!()));
        if t > min_val!() {
            result = result + ((p2 / two!()) * (two!() * phi).sin() / delta).atan() * t.sqrt();
        } else {
            result = result
                + (((p2 / two!()) * (two!() * phi).sin() / delta).atan()
                    * ((one!() / (m - n)).abs()).sqrt()
                    * ((n / (n - one!()).abs()).sqrt()));
        }
        return Ok(result);
    }

    if m == one!() {
        // See http://functions.wolfram.com/08.06.03.0013.01
        result = n.sqrt() * (n.sqrt() * phi.sin()).atanh() - (one!() / phi.cos() + phi.tan()).ln();
        result = result / (n - one!());
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
    let y = one!() - m * t;
    let z = one!();

    let p = if n * t < half!() {
        one!() - n * t
    } else {
        x + vc * t
    };

    let result = sphi * (elliprf(x, y, z)? + n * t * elliprj(x, y, z, p)? / three!());

    Ok(result)
}

#[cfg(test)]
mod tests {
    use crate::compare_test_data;

    use super::*;

    fn ellippiinc_k(inp: &[f64]) -> f64 {
        ellippiinc(inp[1], inp[0], inp[2] * inp[2]).unwrap()
    }

    #[test]
    fn test_ellippi() {
        compare_test_data!(
            "./tests/data/boost/ellippiinc_data.txt",
            ellippiinc_k,
            2.1e-15
        );
    }

    #[test]
    fn test_ellippi_large() {
        compare_test_data!(
            "./tests/data/boost/ellippi3_large_data.txt",
            ellippiinc_k,
            6e-15
        );
    }
}
