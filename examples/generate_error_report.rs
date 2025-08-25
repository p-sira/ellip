/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::*;
use ellip_dev_utils::{env::get_env, get_entry, test_report::generate_error_table};

fn main() {
    let [rust_version, platform, ellip_version] = get_env();
    let template =
        std::fs::read_to_string("examples/error_report_template.md").expect("Cannot read template");

    let env = format!(
        "This report is generated on {} rustc {} using ellip v{} at `f64` precision (ε≈2.22e-16).",
        platform, rust_version, ellip_version
    );
    let legendre_complete = generate_error_table(&[
        get_entry! {"wolfram/ellipk_data", "ellipk", ellipk, 1, 1},
        get_entry! {"wolfram/ellipk_neg", "ellipk (Neg m)", ellipk, 1, 1},
        get_entry! {"wolfram/ellipe_data", "ellipe", ellipe, 1, 1},
        get_entry! {"wolfram/ellipe_neg", "ellipe (Neg m)", ellipe, 1, 1},
        get_entry! {"wolfram/ellippi_data", "ellippi", ellippi, 2, 1},
        get_entry! {"wolfram/ellippi_neg", "ellippi (Neg m)", ellippi, 2, 1},
        get_entry! {"wolfram/ellippi_pv", "ellippi (p.v.)", ellippi, 2, 50},
        get_entry! {"wolfram/ellipd_data", "ellipd", ellipd, 1, 1},
        get_entry! {"wolfram/ellipd_neg", "ellipd (Neg m)", ellipd, 1, 1},
    ]);
    let legendre_incomplete = generate_error_table(&[
        get_entry! {"wolfram/ellipf_data", "ellipf", ellipf, 2, 1},
        get_entry! {"wolfram/ellipf_neg", "ellipf (Neg m)", ellipf, 2, 1},
        get_entry! {"wolfram/ellipeinc_data", "ellipeinc", ellipeinc, 2, 1},
        get_entry! {"wolfram/ellipeinc_neg", "ellipeinc (Neg m)", ellipeinc, 2, 1},
        get_entry! {"wolfram/ellippiinc_data", "ellippiinc", ellippiinc, 3, 1},
        get_entry! {"wolfram/ellippiinc_neg", "ellippiinc (Neg m)", ellippiinc, 3, 1},
        get_entry! {"wolfram/ellippiinc_pv", "ellippiinc (p.v.)", ellippiinc, 3, 1},
        get_entry! {"wolfram/ellippiinc_data", "ellippiinc_bulirsch", ellippiinc_bulirsch, 3, 1},
        get_entry! {"wolfram/ellippiinc_neg", "ellippiinc_bulirsch (Neg m)", ellippiinc_bulirsch, 3, 1},
        get_entry! {"wolfram/ellipdinc_data", "ellipdinc", ellipdinc, 2, 1},
        get_entry! {"wolfram/ellipdinc_neg", "ellipdinc (Neg m)", ellipdinc, 2, 1},
    ]);
    let bulirsch = generate_error_table(&[
        get_entry! {"wolfram/cel_data", "cel", cel, 4, 1},
        get_entry! {"wolfram/cel_pv", "cel (p.v.)", cel, 4, 1},
        get_entry! {"wolfram/cel1_data", "cel1", cel1, 1, 1},
        get_entry! {"wolfram/cel2_data", "cel2", cel2, 3, 1},
        get_entry! {"wolfram/el1_data", "el1", el1, 2, 1},
        get_entry! {"wolfram/el2_data", "el2", el2, 4, 1},
        get_entry! {"wolfram/el3_data", "el3", el3, 3, 50},
        get_entry! {"wolfram/el3_pv", "el3 (p.v.)", el3, 3, 50},
    ]);
    let carlson = generate_error_table(&[
        get_entry! {"wolfram/elliprf_data", "elliprf", elliprf, 3, 1},
        get_entry! {"wolfram/elliprg_data", "elliprg", elliprg, 3, 1},
        get_entry! {"wolfram/elliprj_data", "elliprj", elliprj, 4, 50},
        get_entry! {"wolfram/elliprj_pv", "elliprj (p.v.)", elliprj, 4, 10000000000},
        get_entry! {"boost/elliprj_pv_small", "elliprj (p.v., small*)", elliprj, 4, 50},
        get_entry! {"wolfram/elliprc_data", "elliprc", elliprc, 2, 1},
        get_entry! {"wolfram/elliprc_pv", "elliprc (p.v.)", elliprc, 2, 1},
        get_entry! {"wolfram/elliprd_data", "elliprd", elliprd, 3, 50},
    ]);

    let output = template
        .replace("{{ENV}}", &env)
        .replace("{{LEGENDRE_COMPLETE}}", &legendre_complete)
        .replace("{{LEGENDRE_INCOMPLETE}}", &legendre_incomplete)
        .replace("{{BULIRSCH}}", &bulirsch)
        .replace("{{CARLSON}}", &carlson);

    use std::fs::File;
    use std::io::Write;
    let path = "tests/README.md";
    let mut file = File::create(path).unwrap();
    file.write_all(output.as_bytes()).unwrap();

    println!("Error report generated: {path}");
}
