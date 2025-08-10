/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::*;
use ellip_dev_utils::{env::get_env, get_entry, test_report::generate_error_table};

fn main() {
    let [rust_version, platform, ellip_version] = get_env();

    let lines = [
         "# Testing",
         "This report presents the accuracy of the ellip crate using [**symmetric relative error**](https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/relative_error.html)",
         "metric. Errors are expressed in units of machine epsilon (ε).",
         "The test data spans the domain of each function up to **μ** to avoid approaching the function's limit.",
         "The reference values are computed using [**Wolfram Engine**](https://www.wolfram.com/engine/).",
         "You can find the scripts in the directory [tests/wolfram/](https://github.com/p-sira/ellip/blob/main/tests/wolfram/).",
         &format!("This report is generated on {} rustc {} using ellip v{} at `f64` precision (ε≈2.22e-16).", platform, rust_version, ellip_version),
         "",
         "## Legendre's Complete Elliptic Integrals",
         &generate_error_table(&[
             get_entry!("wolfram/ellipk_data", "ellipk", ellipk, 1, 1),
             get_entry!("wolfram/ellipk_neg", "ellipk (Neg m)", ellipk, 1, 1),
             get_entry!("wolfram/ellipe_data", "ellipe", ellipe, 1, 1),
             get_entry!("wolfram/ellipe_neg", "ellipe (Neg m)", ellipe, 1, 1),
             get_entry!("wolfram/ellippi_data", "ellippi", ellippi, 2, 1),
             get_entry!("wolfram/ellippi_neg", "ellippi (Neg m)", ellippi, 2, 1),
             get_entry!("wolfram/ellippi_pv", "ellippi (p.v.)", ellippi, 2, 50),
         ]),
         "",
         "## Legendre's Incomplete Elliptic Integrals",
         &generate_error_table(&[
             get_entry!("wolfram/ellipf_data", "ellipf", ellipf, 2, 1),
             get_entry!("wolfram/ellipf_neg", "ellipf (Neg m)", ellipf, 2, 1),
             get_entry!("wolfram/ellipeinc_data", "ellipeinc", ellipeinc, 2, 1),
             get_entry!("wolfram/ellipeinc_neg", "ellipeinc (Neg m)", ellipeinc, 2, 1),
             get_entry!("wolfram/ellippiinc_data", "ellippiinc", ellippiinc, 3, 1),
             get_entry!("wolfram/ellippiinc_neg", "ellippiinc (Neg m)", ellippiinc, 3, 1),
             get_entry!("wolfram/ellippiinc_pv", "ellippiinc (p.v.)", ellippiinc, 3, 1),
         ]),
         "",
         "## Bulirsch's Elliptic Integrals",
         "Bulirsh's elliptic integrals are not natively implemented in Wolfram Engine.",
         "Nevertheless, some of the integrals can be converted to their Legendre's counterpart,",
         "which are available on Wolfram Engine. However, for `cel` and `el2`, their values cannot",
         "be mapped simply. Hence, the reference values are generated using the functions",
         "submitted by Jan Mangaldan on [Wolfram Function Repository](https://resources.wolframcloud.com/FunctionRepository/).",
         "As for `cel2`, it is mapped to `cel` with p=1.",
         &generate_error_table(&[
             get_entry!("wolfram/cel_data", "cel", cel, 4, 1),
             get_entry!("wolfram/cel_pv", "cel (p.v.)", cel, 4, 1),
             get_entry!("wolfram/cel1_data", "cel1", cel1, 1, 1),
             get_entry!("wolfram/cel2_data", "cel2", cel2, 3, 1),
             get_entry!("wolfram/el1_data", "el1", el1, 2, 1),
             get_entry!("wolfram/el2_data", "el2", el2, 4, 1),
             get_entry!("wolfram/el3_data", "el3", el3, 3, 50),
             get_entry!("wolfram/el3_pv", "el3 (p.v.)", el3, 3, 50),
         ]),
         "",
         "## Carlson's Symmetric Elliptic Integrals",
         "",
         &generate_error_table(&[
             get_entry!("wolfram/elliprf_data", "elliprf", elliprf, 3, 1),
             get_entry!("wolfram/elliprg_data", "elliprg", elliprg, 3, 1),
             get_entry!("wolfram/elliprj_data", "elliprj", elliprj, 4, 50),
             get_entry!("wolfram/elliprj_pv", "elliprj (p.v.)", elliprj, 4, 10000000000),
             get_entry!("boost/elliprj_pv_small", "elliprj (p.v., small*)", elliprj, 4, 50),
             get_entry!("wolfram/elliprc_data", "elliprc", elliprc, 2, 1),
             get_entry!("wolfram/elliprc_pv", "elliprc (p.v.)", elliprc, 2, 1),
             get_entry!("wolfram/elliprd_data", "elliprd", elliprd, 3, 50),
         ]),
         "",
         "*small: Results compared with Boost Math implementation without promoting double, i.e, computed purely using `f64`.",
         "",
         "Current implementation of `elliprj` is less numerically stable in p.v. cases,",
         "as seen by large errors in the non-small test cases. That said, Ellip's results are consistent",
         "with Boost Math when limited to same precision (See [tests/data/boost/carlson.cpp](https://github.com/p-sira/ellip/blob/main/tests/data/boost/carlson.cpp)).",
         "Since the function is convergent, such errors can be mitigated when Rust's `f128` is released."
     ];

    use std::fs::File;
    use std::io::Write;
    let path = "tests/README.md";
    let mut file = File::create(path).unwrap();

    lines
        .iter()
        .for_each(|line| writeln!(file, "{}", line).expect("Cannot write line"));

    println!("Error report generated: {path}");
}
