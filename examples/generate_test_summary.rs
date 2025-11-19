/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::*;
use ellip_dev_utils::{
    env::{format_clock_speed, get_env},
    get_summary_entry,
    test_report::generate_summary_table,
};

fn main() {
    let env = get_env();
    let template = std::fs::read_to_string("examples/README_template.md")
        .expect("Cannot read README template");

    let summary_section = [
        &format!(
            "Benchmark on {} @{} running {} rustc {} using ellip v{} at `f64` precision (ε=2.2204460492503131e-16).\n",
            env.cpu, format_clock_speed(env.clock_speed), env.platform, env.rust_version, env.ellip_version
        ),
        "### Legendre's Elliptic Integrals",
        &generate_summary_table(&[
            get_summary_entry!("legendre", "ellipk", ellipk, 1),
            get_summary_entry!("legendre", "ellipe", ellipe, 1),
            get_summary_entry!("legendre", "ellipf", ellipf, 2),
            get_summary_entry!("legendre", "ellipeinc", ellipeinc, 2),
            get_summary_entry!("legendre", "ellippi", ellippi, 2),
            get_summary_entry!("legendre", "ellippiinc", ellippiinc, 3),
            get_summary_entry! {"legendre", "ellippiinc_bulirsch", ellippiinc_bulirsch, 3, "ellippiinc"},
            get_summary_entry!("legendre", "ellipd", ellipd, 1),
            get_summary_entry!("legendre", "ellipdinc", ellipdinc, 2),
        ]),
        "",
        "### Bulirsch's Elliptic Integrals",
        &generate_summary_table(&[
            get_summary_entry!("bulirsch", "cel", cel, 4),
            get_summary_entry!("bulirsch", "cel1", cel1, 1),
            get_summary_entry!("bulirsch", "cel2", cel2, 3),
            get_summary_entry!("bulirsch", "el1", el1, 2),
            get_summary_entry!("bulirsch", "el2", el2, 4),
            get_summary_entry!("bulirsch", "el3", el3, 3),
        ]),
        "",
        "### Carlson's Symmetric Integrals",
        &generate_summary_table(&[
            get_summary_entry!("carlson", "elliprf", elliprf, 3),
            get_summary_entry!("carlson", "elliprg", elliprg, 3),
            get_summary_entry!("carlson", "elliprj", elliprj, 4),
            get_summary_entry!("carlson", "elliprc", elliprc, 2),
            get_summary_entry!("carlson", "elliprd", elliprd, 3),
        ]),
        "",
        "### Miscellaneous Functions",
        &generate_summary_table(&[
            get_summary_entry!("misc", "jacobi_zeta", jacobi_zeta, 2),
            get_summary_entry!("misc", "heuman_lambda", heuman_lambda, 2),
        ]),
    ];

    let output = template.replace("{{TESTSUMMARY}}", &summary_section.join("\n"));

    use std::fs::File;
    use std::io::Write;
    let path = "README.md";
    let mut file = File::create(path).unwrap();
    file.write_all(output.as_bytes()).unwrap();

    println!("README generated: {path}");
}
