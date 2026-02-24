/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::*;
use ellip_dev_utils::{
    env::{format_cpu_with_clock_speed, get_env},
    get_entry,
    test_report::generate_error_table,
};

fn main() {
    let env = get_env();
    let template =
        std::fs::read_to_string("examples/error_report_template.md").expect("Cannot read template");

    let env_str = format!(
        "This report is generated on {} running `{} rustc {}` using ellip v{} with `libm` at `f64` precision (ε≈2.22e-16).",
        format_cpu_with_clock_speed(&env.cpu, env.clock_speed), env.platform, env.rust_version, env.ellip_version
    );
    let env_str_f32 = format!(
        "Generated on {} running `{} rustc {}` using ellip v{} with `libm` at `f32` precision (ε≈1.19e-7).",
        format_cpu_with_clock_speed(&env.cpu, env.clock_speed),
        env.platform,
        env.rust_version,
        env.ellip_version
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
        get_entry! {"wolfram/elliprj_pv", "elliprj (p.v.)", elliprj, 4, 1000},
        get_entry! {"wolfram/elliprc_data", "elliprc", elliprc, 2, 1},
        get_entry! {"wolfram/elliprc_pv", "elliprc (p.v.)", elliprc, 2, 1},
        get_entry! {"wolfram/elliprd_data", "elliprd", elliprd, 3, 50},
    ]);
    let misc = generate_error_table(&[
        get_entry! {"wolfram/jacobi_zeta_data", "jacobi_zeta", jacobi_zeta, 2, 1},
        get_entry! {"wolfram/jacobi_zeta_neg", "jacobi_zeta (Neg m)", jacobi_zeta, 2, 1},
        get_entry! {"wolfram/heuman_lambda_data", "heuman_lambda", heuman_lambda, 2, 1},
    ]);

    let output = template
        .replace("{{ENV}}", &env_str)
        .replace("{{LEGENDRE_COMPLETE}}", &legendre_complete[0])
        .replace("{{LEGENDRE_INCOMPLETE}}", &legendre_incomplete[0])
        .replace("{{BULIRSCH}}", &bulirsch[0])
        .replace("{{CARLSON}}", &carlson[0])
        .replace("{{MISC}}", &misc[0])
        .replace("{{ENV_F32}}", &env_str_f32)
        .replace("{{LEGENDRE_COMPLETE_F32}}", &legendre_complete[1])
        .replace("{{LEGENDRE_INCOMPLETE_F32}}", &legendre_incomplete[1])
        .replace("{{BULIRSCH_F32}}", &bulirsch[1])
        .replace("{{CARLSON_F32}}", &carlson[1])
        .replace("{{MISC_F32}}", &misc[1]);

    use std::fs::File;
    use std::io::Write;
    let path = "tests/README.md";
    let mut file = File::create(path).unwrap();
    file.write_all(output.as_bytes()).unwrap();

    println!("Error report generated: {path}");
}
