/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::*;
use ellip_dev_utils::{
    env::get_env, get_accuracy_entry, test_report::generate_accuracy_summary_table,
};

// Generate accuracy summary in the JOSS article style
fn main() {
    let env = get_env();

    let output = [
        "# Accuracy summary",
        "",
        "Accuracy summary in the JOSS article style.",
        &format!(
            "Generated on {} rustc {} using ellip v{} at `f64` precision (ε=2.2204460492503131e-16).\n",
            env.platform, env.rust_version, env.ellip_version
        ),
        &generate_accuracy_summary_table(&[
            get_accuracy_entry!("legendre", "ellipk", ellipk, 1),
            get_accuracy_entry!("legendre", "ellipe", ellipe, 1),
            get_accuracy_entry!("legendre", "ellippi", ellippi, 2),
            get_accuracy_entry!("legendre", "ellipd", ellipd, 1),
            get_accuracy_entry!("legendre", "ellipf", ellipf, 2),
            get_accuracy_entry!("legendre", "ellipeinc", ellipeinc, 2),
            get_accuracy_entry!("legendre", "ellippiinc", ellippiinc, 3),
            get_accuracy_entry! {"legendre", "ellippiinc_bulirsch", ellippiinc_bulirsch, 3, "ellippiinc"},
            get_accuracy_entry!("legendre", "ellipdinc", ellipdinc, 2),
            get_accuracy_entry!("bulirsch", "cel", cel, 4),
            get_accuracy_entry!("bulirsch", "cel1", cel1, 1),
            get_accuracy_entry!("bulirsch", "cel2", cel2, 3),
            get_accuracy_entry!("bulirsch", "el1", el1, 2),
            get_accuracy_entry!("bulirsch", "el2", el2, 4),
            get_accuracy_entry!("bulirsch", "el3", el3, 3),
            get_accuracy_entry!("carlson", "elliprf", elliprf, 3),
            get_accuracy_entry!("carlson", "elliprg", elliprg, 3),
            get_accuracy_entry!("carlson", "elliprj", elliprj, 4),
            get_accuracy_entry!("carlson", "elliprc", elliprc, 2),
            get_accuracy_entry!("carlson", "elliprd", elliprd, 3),
            get_accuracy_entry!("misc", "jacobi_zeta", jacobi_zeta, 2),
            get_accuracy_entry!("misc", "heuman_lambda", heuman_lambda, 2),
        ]),
    ].join("\n");

    use std::fs::File;
    use std::io::Write;
    let path = "tests/accuracy_summary.md";
    let mut file = File::create(path).unwrap();
    file.write_all(output.as_bytes()).unwrap();

    println!("Accuracy summary generated: {path}");
}
