/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

#[macro_export]
macro_rules! func_wrapper {
    ($func:expr, 1) => {
        fn wrapped_func(args: &Vec<f64>) -> f64 {
            $func(args[0]).unwrap()
        }
    };
    ($func:expr, 2) => {
        fn wrapped_func(args: &Vec<f64>) -> f64 {
            $func(args[0], args[1]).unwrap()
        }
    };
    ($func:expr, 3) => {
        fn wrapped_func(args: &Vec<f64>) -> f64 {
            let ans = $func(args[0], args[1], args[2]).unwrap_or_else(|e| {
                println!("{}: {}, {}, {}", e, args[0], args[1], args[2]);
                panic!()
            });
            if ans.is_nan() {
                println!("{}, {}, {}", args[0], args[1], args[2]);
                panic!("NAN")
            }
            ans
        }
    };
    ($func:expr, 4) => {
        fn wrapped_func(args: &Vec<f64>) -> f64 {
            $func(args[0], args[1], args[2], args[3]).unwrap()
        }
    };
    ($func:expr, $_:expr) => {
        fn wrapped_func(_args: &Vec<f64>) -> f64 {
            panic!("Unsupported number of arguments")
        }
    };
}
