/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

//! # Ellip Rayon
//! Parallelized elliptic integral computation for Rust based on [Ellip](https://github.com/p-sira/ellip).
//!
//! ## Installation
//! ```shell
//! cargo add ellip-rayon
//! ```
//!
//! ## Machine-specific Threshold
//! Ellip Rayon employs parallelization if the length of the arguments exceesds certain thresholds. These thresholds depend on the core count, cache size, and architecture. For the most efficiency, these thresholds should be tuned on the target machine.
//!
//! 1. Generating benchmark data
//! 
//! ```shell
//! cargo bench
//! ```
//! The process should take about 30-40 minutes to complete, and the file `benches/par_threshold.md` will be created. This reports the threshold for each function.
//!
//! 2. Using the generated threshold
//! 
//! ```shell
//! cargo run --example generate_threshold_code
//! ```
//! This script automatically replaces the thresholds in the source code.
//!
//! 3. Adding locally compiled library
//! 
//! From your working directory, run
//! ```shell
//! cargo add --path path/to/your/ellip-rayon
//! ```

use ellip::*;
use itertools::izip;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

macro_rules! par_zip {
    ($a:expr) => {
        $a.par_iter()
    };
    ($a:expr, $b:expr) => {
        $a.par_iter().zip($b.par_iter())
    };
    ($a:expr, $b:expr, $c:expr) => {
        $a.par_iter()
            .zip($b.par_iter())
            .zip($c.par_iter())
            .map(|((a, b), c)| (a, b, c))
    };
    ($a:expr, $b:expr, $c:expr, $d:expr) => {
        $a.par_iter()
            .zip($b.par_iter())
            .zip($c.par_iter())
            .zip($d.par_iter())
            .map(|(((a, b), c), d)| (a, b, c, d))
    };
}

macro_rules! impl_par {
    (@inner, $fn:ident, 2) => {
        |(&a, &b)| ellip::$fn(a, b)
    };
    (@inner, $fn:ident, 3) => {
        |(&a, &b, &c)| ellip::$fn(a, b, c)
    };
    (@inner, $fn:ident, 4) => {
        |(&a, &b, &c, &d)| ellip::$fn(a, b, c, d)
    };
    ($fn:ident, [$arg:ident], 1) => {
        #[doc=concat!["Computes [", stringify!($fn), "](ellip::", stringify!($fn), ") in parallel."]]
        pub fn $fn($arg: &[f64]) -> Result<Vec<f64>, StrErr> {
            $arg.iter().map(|&a| ellip::$fn(a)).collect()
        }
    };
    ($fn:ident, [$arg:ident], 1, $threshold:expr) => {
        #[doc=concat!["Computes [", stringify!($fn), "](ellip::", stringify!($fn), ") in parallel."]]
        pub fn $fn($arg: &[f64]) -> Result<Vec<f64>, StrErr> {
            if $arg.len() < $threshold {
                $arg.iter().map(|&a| ellip::$fn(a)).collect()
            } else {
                $arg.par_iter().map(|&a| ellip::$fn(a)).collect()
            }
        }
    };
    ($fn:ident, [$first:ident, $($args:ident),*], $n_arg:tt) => {
        #[doc=concat!["Computes [", stringify!($fn), "](ellip::", stringify!($fn), ") in parallel."]]
        pub fn $fn($first: &[f64], $($args: &[f64],)*) -> Result<Vec<f64>, StrErr> {
            $(
                if $first.len() != $args.len() {
                    return Err(concat![stringify!($fn), ": All arguments must have the same length."]);
                }
            )*
            izip!($first, $($args),*).map(impl_par!(@inner, $fn, $n_arg)).collect()
        }
    };
    ($fn:ident, [$first:ident, $($args:ident),*], $n_arg:tt, $threshold:expr) => {
        #[doc=concat!["Computes [", stringify!($fn), "](ellip::", stringify!($fn), ") in parallel."]]
        pub fn $fn($first: &[f64], $($args: &[f64],)*) -> Result<Vec<f64>, StrErr> {
            $(
                if $first.len() != $args.len() {
                    return Err(concat![stringify!($fn), ": All arguments must have the same length."]);
                }
            )*
            if $first.len() < $threshold {
                izip!($first, $($args),*).map(impl_par!(@inner, $fn, $n_arg)).collect()
            } else {
                par_zip!($first, $($args),*).map(impl_par!(@inner, $fn, $n_arg)).collect()
            }
        }
    };
}

// {BEGIN_IMPL_PAR}
// Generated threshold values from benchmark results
impl_par!(ellipk, [m], 1, 1000);
impl_par!(ellipe, [m], 1, 1500);
impl_par!(ellipf, [phi, m], 2, 400);
impl_par!(ellipeinc, [phi, m], 2, 300);
impl_par!(ellippi, [n, m], 2, 200);
impl_par!(ellippiinc, [phi, n, m], 3, 200);
impl_par!(ellippiinc_bulirsch, [phi, n, m], 3, 300);
impl_par!(ellipd, [m], 1, 600);
impl_par!(ellipdinc, [phi, m], 2, 500);
impl_par!(cel, [kc, p, a, b], 4, 600);
impl_par!(cel1, [kc], 1, 1500);
impl_par!(cel2, [kc, a, b], 3, 700);
impl_par!(el1, [x, kc], 2, 600);
impl_par!(el2, [x, kc, a, b], 4, 600);
impl_par!(el3, [x, kc, p], 3, 500);
impl_par!(elliprf, [x, y, z], 3, 600);
impl_par!(elliprg, [x, y, z], 3, 500);
impl_par!(elliprj, [x, y, z, p], 4, 300);
impl_par!(elliprc, [x, y], 2, 800);
impl_par!(elliprd, [x, y, z], 3, 500);
impl_par!(jacobi_zeta, [phi, m], 2, 200);
impl_par!(heuman_lambda, [phi, m], 2, 200);

// {END_IMPL_PAR}
