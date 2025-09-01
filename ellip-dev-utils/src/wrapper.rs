/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

#[macro_export]
macro_rules! func_wrapper {
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
