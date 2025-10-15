/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

pub struct Stats {
    pub mean: f64,
    pub median: f64,
    pub variance: f64,
    pub max: f64,
    pub p99: f64,
    pub n: usize,
}

impl Stats {
    pub fn nan() -> Self {
        return Stats {
            mean: f64::NAN,
            median: f64::NAN,
            variance: f64::NAN,
            max: f64::NAN,
            p99: f64::NAN,
            n: 0,
        };
    }

    pub fn from_vec(v: &Vec<f64>) -> Self {
        let mut valids: Vec<f64> = v
            .iter()
            .filter(|x| !x.is_nan() && x.is_finite())
            .cloned()
            .collect();

        if valids.is_empty() {
            return Self::nan();
        }

        valids.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let sum: f64 = valids.iter().sum();
        let n = valids.len();
        let mean = sum / n as f64;

        // Calculate median
        let median = if n % 2 == 0 {
            (valids[n / 2 - 1] + valids[n / 2]) / 2.0
        } else {
            valids[n / 2]
        };

        // Calculate P99 error
        let p99_pos = (n - 1) as f64 * 0.99;
        let p99_pos_low = p99_pos.floor() as usize;
        let p99_frac = p99_pos - p99_pos_low as f64;

        let p99 = if p99_pos_low + 1 < n {
            valids[p99_pos_low] * (1.0 - p99_frac) + valids[p99_pos_low + 1] * p99_frac
        } else {
            valids[n - 1]
        };

        // Calculate variance
        let variance = valids
            .iter()
            .map(|x| {
                let diff = mean - x;
                diff * diff
            })
            .sum::<f64>()
            / n as f64;

        let max = *valids
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap_or(&f64::NAN);

        Stats {
            mean,
            median,
            variance,
            max,
            p99,
            n,
        }
    }
}
