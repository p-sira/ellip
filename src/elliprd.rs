/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 * This code is modified from Russell Lab to achieve precision close to SciPy.
 */

// Original header from Russell Lab
// Computes elliptic integral of the second kind using Carlson's formula
//
// Computes Rd(x,y,z) where x,y must be non-negative and at most one can be zero. z must be positive.
//
// # References:
//
// * Press WH, Teukolsky SA, Vetterling WT, Flannery BP (2007) Numerical Recipes: The Art of
//   Scientific Computing. Third Edition. Cambridge University Press. 1235p.
// * Carlson BC (1977) Elliptic Integrals of the First Kind, SIAM Journal on Mathematical
//   Analysis, vol. 8, pp. 231-242.

/// Compute degenerate symmetric elliptic integral of the third kind.
/// ```text
///                     ∞                                        
///                 3  ⌠                   dt                   
/// RD(x, y, z)  =  ─  ⎮ ───────────────────────────────────────
///                 2  ⎮             ___________________________
///                    ⌡ (t + z) ⋅ ╲╱(t + x) ⋅ (t + y) ⋅ (t + z)
///                   0                                        
/// where x ≥ 0, y ≥ 0, and at most one can be zero. z > 0.
/// ```
///
pub fn elliprd(x: f64, y: f64, z: f64) -> Result<f64, &'static str> {
    let tiny = 2.0 * f64::MAX.powf(-2.0 / 3.0);
    let big = 0.1 * RD_ERR_TOL * f64::MIN_POSITIVE.powf(-2.0 / 3.0);
    if x.min(y) < 0.0 || (x + y).min(z) < tiny || x.max(y).max(z) > big {
        return Err("elliprf: x and y must be non-negative, and at most one can be zero. z must be positive");
    }
    let mut xt = x;
    let mut yt = y;
    let mut zt = z;
    let mut ave: f64 = 0.0;
    let mut dx: f64 = 0.0;
    let mut dy: f64 = 0.0;
    let mut dz: f64 = 0.0;
    let mut sum = 0.0;
    let mut fac = 1.0;
    let mut it = 0;
    for _ in 0..N_MAX_ITERATIONS {
        let sqx = xt.sqrt();
        let sqy = yt.sqrt();
        let sqz = zt.sqrt();
        let lam = sqx * (sqy + sqz) + sqy * sqz;
        sum += fac / (sqz * (zt + lam));
        fac *= 0.25;
        xt = 0.25 * (xt + lam);
        yt = 0.25 * (yt + lam);
        zt = 0.25 * (zt + lam);
        ave = 0.2 * (xt + yt + 3.0 * zt);
        dx = (ave - xt) / ave;
        dy = (ave - yt) / ave;
        dz = (ave - zt) / ave;
        if dx.abs().max(dy.abs()).max(dz.abs()) < RD_ERR_TOL {
            break;
        }
        it += 1;
    }
    if it == N_MAX_ITERATIONS {
        return Err("elliprd: Fail to converge.");
    }
    let ea = dx * dy;
    let eb = dz * dz;
    let ec = ea - eb;
    let ed = ea - 6.0 * eb;
    let ee = ed + ec + ec;
    let ans = 3.0 * sum
        + fac
            * (1.0
                + ed * (-RD_C1 + RD_C5 * ed - RD_C6 * dz * ee)
                + dz * (RD_C2 * ee + dz * (-RD_C3 * ec + dz * RD_C4 * ea)))
            / (ave * ave.sqrt());
    Ok(ans)
}

const N_MAX_ITERATIONS: usize = 500; // Modified from Russell's 11
const RD_ERR_TOL: f64 = 0.0015;
const RD_C1: f64 = 3.0 / 14.0;
const RD_C2: f64 = 1.0 / 6.0;
const RD_C3: f64 = 9.0 / 22.0;
const RD_C4: f64 = 3.0 / 26.0;
const RD_C5: f64 = 0.25 * RD_C3;
const RD_C6: f64 = 1.5 * RD_C4;

#[cfg(test)]
mod test {
    use super::*;
    use crate::compare_test_data;

    fn _elliprd(inp: &[f64]) -> f64 {
        elliprd(inp[0], inp[1], inp[2]).unwrap()
    }

    #[test]
    fn test_elliprd() {
        compare_test_data!("./tests/data/boost/ellint_rd_data.txt", _elliprd, 6e-16);
    }
}
