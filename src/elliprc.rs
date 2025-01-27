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

use crate::constants::{BIG, TINY};

/// Compute degenerate symmetric elliptic integral of RF.
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
pub fn elliprc(x: f64, y: f64) -> Result<f64, &'static str> {
    let comp1 = 2.236 / TINY.sqrt();
    let comp2 = (TINY * BIG).powi(2) / 25.0;
    if x < 0.0
        || y == 0.0
        || (x + y.abs()) < TINY
        || (x + y.abs()) > BIG
        || (y < -comp1 && x > 0.0 && x < comp2)
    {
        return Err("elliprc: input must satisfy: x ≥ 0, y ≠ 0.");
    }
    let (mut xt, mut yt, w) = if y > 0.0 {
        (x, y, 1.0)
    } else {
        (x - y, -y, x.sqrt() / (x - y).sqrt())
    };
    let mut ave: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut it = 0;
    for _ in 0..N_MAX_ITERATIONS {
        let lam = 2.0 * xt.sqrt() * yt.sqrt() + yt;
        xt = 0.25 * (xt + lam);
        yt = 0.25 * (yt + lam);
        ave = (xt + yt + yt) / 3.0;
        s = (yt - ave) / ave;
        if s.abs() < RC_ERR_TOL {
            break;
        }
        it += 1;
    }
    if it == N_MAX_ITERATIONS {
        return Err("elliprc: Fail to converge.");
    }
    let ans = w * (1.0 + s * s * (RC_C1 + s * (RC_C2 + s * (RC_C3 + s * RC_C4)))) / ave.sqrt();
    Ok(ans)
}

const N_MAX_ITERATIONS: usize = 500; // Modified from Russell's 11
const RC_ERR_TOL: f64 = 0.0006; // Modified from Russell's 0.0012
const RC_C1: f64 = 0.3;
const RC_C2: f64 = 1.0 / 7.0;
const RC_C3: f64 = 0.375;
const RC_C4: f64 = 9.0 / 22.0;

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
