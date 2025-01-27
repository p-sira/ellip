/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Russell Lab to achieve precision on par with SciPy.
 */

// Original header from Russell Lab
// Computes elliptic integral of the first kind using Carlson's formula
//
// Computes Rf(x,y,z) where x,y,z must be non-negative and at most one can be zero.
//
// # References:
//
// * Press WH, Teukolsky SA, Vetterling WT, Flannery BP (2007) Numerical Recipes: The Art of
//   Scientific Computing. Third Edition. Cambridge University Press. 1235p.
// * Carlson BC (1977) Elliptic Integrals of the First Kind, SIAM Journal on Mathematical
//   Analysis, vol. 8, pp. 231-242.

use crate::constants::{BIG, TINY};

/// Compute symmetric elliptic integral of the first kind.
/// ```text
///                     ∞                              
///                 1  ⌠              dt              
/// RF(x, y, z)  =  ─  ⎮ ─────────────────────────────
///                 2  ⎮   ___________________________
///                    ⌡ ╲╱(t + x) ⋅ (t + y) ⋅ (t + z)
///                   0                              
/// where x ≥ 0, y ≥ 0, z ≥ 0, and at most one can be zero.
/// ```
///
pub fn elliprf(x: f64, y: f64, z: f64) -> Result<f64, &'static str> {
    if x.min(y).min(z) < 0.0 || (y + z).min(x + y).min(x + z) < TINY || x.max(y).max(z) > BIG {
        return Err("elliprf: x, y, and z must be non-negative, and at most one can be zero.");
    }
    let mut xt = x;
    let mut yt = y;
    let mut zt = z;
    let mut ave: f64 = 0.0;
    let mut dx: f64 = 0.0;
    let mut dy: f64 = 0.0;
    let mut dz: f64 = 0.0;
    let mut it = 0;
    for _ in 0..N_MAX_ITERATIONS {
        let sqx = xt.sqrt();
        let sqy = yt.sqrt();
        let sqz = zt.sqrt();
        let lam = sqx * (sqy + sqz) + sqy * sqz;
        xt = 0.25 * (xt + lam);
        yt = 0.25 * (yt + lam);
        zt = 0.25 * (zt + lam);
        ave = (xt + yt + zt) / 3.0;
        dx = (ave - xt) / ave;
        dy = (ave - yt) / ave;
        dz = (ave - zt) / ave;
        if dx.abs().max(dy.abs()).max(dz.abs()) < RF_ERR_TOL {
            break;
        }
        it += 1;
    }
    if it == N_MAX_ITERATIONS {
        return Err("elliprf: Fail to converge.");
    }
    let e2 = dx * dy - dz * dz;
    let e3 = dx * dy * dz;
    let ans = (1.0 + (RF_C1 * e2 - RF_C2 - RF_C3 * e3) * e2 + RF_C4 * e3) / ave.sqrt();
    Ok(ans)
}

const N_MAX_ITERATIONS: usize = 500; // Modified from Russell's 11
const RF_ERR_TOL: f64 = 0.0001; // Modified from Russell's 0.0025
const RF_C1: f64 = 1.0 / 24.0;
const RF_C2: f64 = 0.1;
const RF_C3: f64 = 3.0 / 44.0;
const RF_C4: f64 = 1.0 / 14.0;

#[cfg(test)]
mod test {
    use super::*;
    use crate::compare_test_data;

    fn _elliprf(inp: &[f64]) -> f64 {
        elliprf(inp[0], inp[1], inp[2]).unwrap()
    }

    #[test]
    fn test_elliprf() {
        compare_test_data!("./tests/data/boost/ellint_rf_data.txt", _elliprf, 5e-16);
    }
}
