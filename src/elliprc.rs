/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Russell Lab to achieve precision on par with SciPy.
 */

// Original header from Russell Lab
// Computes the degenerate elliptic integral using Carlson's formula
//
// Computes Rc(x,y) where x must be non-negative and y must be nonzero.
// If y < 0, the Cauchy principal value is returned.
//
// # References:
//
// * Press WH, Teukolsky SA, Vetterling WT, Flannery BP (2007) Numerical Recipes: The Art of
//   Scientific Computing. Third Edition. Cambridge University Press. 1235p.

use crate::{BIG, TINY};
use num_traits::Float;

/// Compute [degenerate symmetric elliptic integral of RF](https://dlmf.nist.gov/19.16.E6).
/// ```text
///                  ∞                  
///              1  ⌠        dt        
/// RC(x, y)  =  ─  ⎮ ─────────────────
///              2  ⎮             _____
///                 ⌡ (t + y) ⋅ ╲╱t + x
///                0                  
/// where x ≥ 0, y ≠ 0
/// ```
///
pub fn elliprc<T: Float>(x: T, y: T) -> Result<T, &'static str> {
    let tiny = T::from(TINY).unwrap();
    let big = T::from(BIG).unwrap();
    let comp1 = T::from(2.236).unwrap() / tiny.sqrt();
    let comp2 = (tiny * big).powi(2) / T::from(25.0).unwrap();

    if x < T::zero()
        || y == T::zero()
        || (x + y.abs()) < tiny
        || (x + y.abs()) > big
        || (y < -comp1 && x > T::zero() && x < comp2)
    {
        return Err("elliprc: input must satisfy: x ≥ 0, y ≠ 0.");
    }

    let (mut xt, mut yt, w) = if y > T::zero() {
        (x, y, T::one())
    } else {
        (x - y, -y, x.sqrt() / (x - y).sqrt())
    };

    let n_max_iterations = 500;
    let rc_err_tol = T::from(0.0006).unwrap();
    let rc_c1 = T::from(0.3).unwrap();
    let rc_c2 = T::from(1.0 / 7.0).unwrap();
    let rc_c3 = T::from(0.375).unwrap();
    let rc_c4 = T::from(9.0 / 22.0).unwrap();

    let mut ave: T = T::zero();
    let mut s: T = T::zero();
    let mut it = 0;

    for _ in 0..n_max_iterations {
        let lam = T::from(2.0).unwrap() * xt.sqrt() * yt.sqrt() + yt;
        xt = T::from(0.25).unwrap() * (xt + lam);
        yt = T::from(0.25).unwrap() * (yt + lam);
        ave = (xt + yt + yt) / T::from(3.0).unwrap();
        s = (yt - ave) / ave;

        if s.abs() < rc_err_tol {
            break;
        }
        it += 1;
    }

    if it == n_max_iterations {
        return Err("elliprc: Fail to converge.");
    }

    let ans = w * (T::one() + s * s * (rc_c1 + s * (rc_c2 + s * (rc_c3 + s * rc_c4)))) / ave.sqrt();
    Ok(ans)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::compare_test_data;

    fn _elliprc(inp: &[f64]) -> f64 {
        elliprc(inp[0], inp[1]).unwrap()
    }

    #[test]
    fn test_elliprc() {
        compare_test_data!("./tests/data/boost/ellint_rc_data.txt", _elliprc, 5e-16);
    }
}
