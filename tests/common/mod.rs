/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

#[macro_export]
macro_rules! compare_test_data {
    ($func: expr, $cases :expr, $t: ident, $rtol:expr, $atol:expr) => {
        let mut test_fail = 0;
        let mut worst_line = 0;
        let mut worst_params: Vec<$t> = Vec::new();
        let mut worst_errors: [$t;5] = [0.0; 5];

        for (i, case) in $cases.iter().enumerate() {
            let inp = &case.inputs;
            let exp = case.expected;
            let result = $func(&inp);

            if exp.is_nan() && result.is_nan() {
                continue;
            }

            if result.is_nan() {
                eprintln! ("Test failed on case {}: input = {:?}, expected = {:?}, got = NAN",
                    i + 1,
                    inp,
                    exp,
                );
                test_fail += 1;
            }

            if !result.is_finite() {
                if result != exp {
                    eprintln! ("Test failed on case {}: input = {:?}, expected = {:?}, got = {:?}",
                        i + 1,
                        inp,
                        exp,
                        result);
                    test_fail += 1;
                }
            }

            let error = (result - exp).abs();
            let tol = $t::from($atol) + $t::from($rtol) * exp.abs();
            if error > tol {
                let rel = error / exp.abs();
                test_fail += 1;
                if rel > worst_errors[1] {
                    worst_line = i + 1;
                    worst_errors = [error, rel, tol, exp, result];
                    worst_params = inp.clone();
                }
                eprintln!(
                    "Test failed on case {}: input = {:?}, expected = {:?}, got = {:?} \n error = {:?}, rel = {:?}, abs = {:?}, rtol = {:?}, atol = {:?}\n",
                    i + 1,
                    inp,
                    exp,
                    result,
                    error,
                    rel,
                    tol,
                    $rtol,
                    $atol,
                );
            }
        }

        if test_fail > 0 {
            panic!(
                "Failed {}/{} cases. Worst on case {}: input = {:?}, expected = {:?}, got = {:?} \n error = {:?}, rel = {:?}, abs = {:?}, rtol = {:?}, atol = {:?}",
                test_fail,
                $cases.len(),
                worst_line,
                worst_params,
                worst_errors[3],
                worst_errors[4],
                worst_errors[0],
                worst_errors[1],
                worst_errors[2],
                $rtol,
                $atol,
            );
        }
    };
}

#[macro_export]
macro_rules! compare_test_data_boost {
    ($filename:expr, $func:expr, $rtol:expr) => {
        compare_test_data_boost!($filename, $func, f64, $rtol, 0.0)
    };
    ($filename:expr, $func:expr, $rtol:expr, atol: $atol:expr) => {
        compare_test_data_boost!($filename, $func, f64, $rtol, $atol)
    };
    ($filename:expr, $func:expr, $n_args:tt, $rtol:expr) => {{
        use ellip_dev_utils::func_wrapper;
        func_wrapper!($func, k, $n_args);
        compare_test_data_boost!($filename, wrapped_func, f64, $rtol, 0.0)
    }};
    ($filename:expr, $func:expr, $t:ident, $rtol:expr, $atol:expr) => {{
        use ellip_dev_utils::parser;
        use std::path::Path;

        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/boost")
            .join($filename);
        let cases = parser::read_boost_data(&path.to_str().unwrap()).unwrap_or_else(|err| {
            panic!("Failed to read reference dataset {}: {err}", path.display())
        });
        $crate::compare_test_data!($func, cases, $t, $rtol, $atol);
    }};
}

#[macro_export]
macro_rules! compare_test_data_wolfram {
    ($filename:expr, $func:expr, $n_args:tt, $rtol:expr, atol: $atol:expr) => {{
        use ellip_dev_utils::func_wrapper;
        func_wrapper!($func, $n_args);
        compare_test_data_wolfram!(
            "./tests/data/wolfram",
            $filename,
            wrapped_func,
            f64,
            $rtol,
            $atol
        )
    }};
    ($filename:expr, $func:expr, $n_args:tt, $rtol:expr) => {
        compare_test_data_wolfram!($filename, $func, $n_args, $rtol, atol: 0.0)
    };
    ($path:expr, $filename:expr, $func:expr, $n_args:tt, $rtol:expr) => {{
        use ellip_dev_utils::func_wrapper;
        func_wrapper!($func, $n_args);
        compare_test_data_wolfram!($path, $filename, wrapped_func, f64, $rtol, 0.0)
    }};
    ($path:expr, $filename:expr, $func:expr, $t:ident, $rtol:expr, $atol:expr) => {{
        use ellip_dev_utils::parser;
        use std::path::Path;

        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join($path)
            .join($filename);
        let cases = parser::read_wolfram_data(&path.to_str().unwrap()).unwrap_or_else(|err| {
            panic!("Failed to read reference dataset {}: {err}", path.display())
        });
        $crate::compare_test_data!($func, cases, $t, $rtol, $atol);
    }};
}

pub fn linspace(start: f64, end: f64, num: usize) -> Vec<f64> {
    if num < 2 {
        return vec![start];
    }

    let step = (end - start) / (num - 1) as f64;
    (0..num).map(|i| start + step * i as f64).collect()
}

mod dataset_regression_tests {
    // Regression for audit finding A16: https://github.com/p-sira/ellip/pull/130
    fn identity(input: &[f64]) -> f64 {
        input[0]
    }

    // Regression for audit finding A16: https://github.com/p-sira/ellip/pull/130
    #[test]
    #[should_panic(expected = "missing-audit-a16.txt")]
    fn missing_boost_dataset_fails() {
        crate::compare_test_data_boost!("missing-audit-a16.txt", identity, 1e-15);
    }

    // Regression for audit finding A16: https://github.com/p-sira/ellip/pull/130
    #[test]
    #[should_panic(expected = "missing-audit-a16.csv")]
    fn missing_wolfram_dataset_fails() {
        crate::compare_test_data_wolfram!(
            "tests/data/wolfram",
            "missing-audit-a16.csv",
            identity,
            f64,
            1e-15,
            0.0
        );
    }

    // Regression for audit finding A16: https://github.com/p-sira/ellip/pull/130
    #[test]
    fn datasets_load_outside_checkout() {
        for test in [
            "accuracy::elliprc::test_elliprc",
            "accuracy::ellipk::test_ellipk_wolfram",
            "accuracy::elliprj::test_elliprj_pv_rounding_f32",
        ] {
            let output = std::process::Command::new(std::env::current_exe().unwrap())
                .args(["--exact", test])
                .current_dir(std::env::temp_dir())
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "{test}: {}\n{}",
                String::from_utf8_lossy(&output.stdout),
                String::from_utf8_lossy(&output.stderr)
            );
            let stdout = String::from_utf8(output.stdout).unwrap();
            assert!(stdout.contains("1 passed"), "{test}: {stdout}");
        }
    }
}
