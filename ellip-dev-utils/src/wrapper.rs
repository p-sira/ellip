/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

#[macro_export]
macro_rules! func_wrapper {
    ($func:expr, $n_args: tt) => {
        func_wrapper!($func, f64, $n_args)
    };
    ($func:expr, k,  1) => {
        fn wrapped_func(args: &Vec<f64>) -> f64 {
            $func(args[0] * args[0]).unwrap()
        }
    };
    ($func:expr, k, 2) => {
        fn wrapped_func(args: &Vec<f64>) -> f64 {
            $func(args[0], args[1] * args[1]).unwrap()
        }
    };
    ($func:expr, k, 3) => {
        fn wrapped_func(args: &Vec<f64>) -> f64 {
            $func(args[0], args[1], args[2] * args[2]).unwrap()
        }
    };
    ($func:expr, $t: ty, 1) => {
        fn wrapped_func(args: &Vec<$t>) -> $t {
            $func(args[0]).unwrap()
        }
    };
    ($func:expr, $t: ty, 2) => {
        fn wrapped_func(args: &Vec<$t>) -> $t {
            $func(args[0], args[1]).unwrap()
        }
    };
    ($func:expr, $t: ty, 3) => {
        fn wrapped_func(args: &Vec<$t>) -> $t {
            $func(args[0], args[1], args[2]).unwrap()
        }
    };
    ($func:expr, $t: ty, 4) => {
        fn wrapped_func(args: &Vec<$t>) -> $t {
            $func(args[0], args[1], args[2], args[3]).unwrap()
        }
    };
    ($func:expr, $t: ty, $_:expr) => {
        fn wrapped_func(_args: &Vec<$t>) -> $t {
            panic!("Unsupported number of arguments")
        }
    };
}
