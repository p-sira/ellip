/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use std::fs;
use std::path::Path;

use regex::Regex;

#[derive(Debug, Clone)]
struct BenchmarkResult {
    function_name: String,
    threshold: usize,
    exit_condition: String,
}

fn parse_markdown_table(content: &str) -> Vec<BenchmarkResult> {
    let mut results = Vec::new();
    let lines: Vec<&str> = content.lines().collect();

    // Skip header lines (first 2 lines are | Function | ... | and |-----|...)
    for line in lines.iter().skip(2) {
        if line.trim().is_empty() || !line.starts_with('|') {
            continue;
        }

        let parts: Vec<&str> = line.split('|').collect();
        if parts.len() >= 6 {
            let function_name = parts[1].trim().to_string();
            let threshold_str = parts[2].trim();
            let exit_condition = parts[6].trim().to_string();

            if let Ok(threshold) = threshold_str.parse::<usize>() {
                results.push(BenchmarkResult {
                    function_name,
                    threshold,
                    exit_condition,
                });
            }
        }
    }

    results
}

fn generate_code(func_name: &str, func_data: &str, bench_result: BenchmarkResult) -> String {
    if bench_result.exit_condition == "Min step size" || bench_result.exit_condition == "Converged"
    {
        format!(
            "impl_par!({}, {}, {});\n",
            func_name, func_data, bench_result.threshold
        )
    } else {
        format!("impl_par!({}, {});\n", func_name, func_data)
    }
}

fn generate_lib_rs_code(results: &[BenchmarkResult]) -> String {
    let mut code = String::new();

    code.push_str("// Generated threshold values from benchmark results\n");

    let functions = [
        // Legendre's Integrals
        ("ellipk", "[m], 1"),
        ("ellipe", "[m], 1"),
        ("ellipf", "[phi, m], 2"),
        ("ellipeinc", "[phi, m], 2"),
        ("ellippi", "[n, m], 2"),
        ("ellippiinc", "[phi, n, m], 3"),
        ("ellippiinc_bulirsch", "[phi, n, m], 3"),
        ("ellipd", "[m], 1"),
        ("ellipdinc", "[phi, m], 2"),
        // Bulirsch's Integrals
        ("cel", "[kc, p, a, b], 4"),
        ("cel1", "[kc], 1"),
        ("cel2", "[kc, a, b], 3"),
        ("el1", "[x, kc], 2"),
        ("el2", "[x, kc, a, b], 4"),
        ("el3", "[x, kc, p], 3"),
        // Carlson's Integrals
        ("elliprf", "[x, y, z], 3"),
        ("elliprg", "[x, y, z], 3"),
        ("elliprj", "[x, y, z, p], 4"),
        ("elliprc", "[x, y], 2"),
        ("elliprd", "[x, y, z], 3"),
        // Miscellaneous Functions
        ("jacobi_zeta", "[phi, m], 2"),
        ("heuman_lambda", "[phi, m], 2"),
    ];

    for func in &functions {
        if let Some(result) = results.iter().find(|r| r.function_name == *func.0) {
            code.push_str(&generate_code(func.0, func.1, result.clone()));
        } else {
            panic!("Benchmark not found for function {}", func.0)
        }
    }

    code
}

fn main() {
    let md_file_path = Path::new("benches/par_threshold.md");

    if !md_file_path.exists() {
        eprintln!("Error: {} not found", md_file_path.display());
        panic!("Please run the benchmark first to generate the results.");
    }

    let content = match fs::read_to_string(md_file_path) {
        Ok(content) => content,
        Err(e) => {
            panic!("Error reading {}: {}", md_file_path.display(), e);
        }
    };

    let results = parse_markdown_table(&content);

    if results.is_empty() {
        panic!("No benchmark results found in {}", md_file_path.display());
    }

    let generated_code = generate_lib_rs_code(&results);
    println!("{generated_code}");

    let lib_path = Path::new("src/lib.rs");

    if !lib_path.exists() {
        eprintln!("Error: {} not found", lib_path.display());
    }

    let content = fs::read_to_string(lib_path).expect("Cannot read lib.rs");

    let re = Regex::new(r"(?s)// \{\{BEGIN_IMPL_PAR\}\}.*?// \{\{END_IMPL_PAR\}\}")
        .expect("Invalid regex");
    let output = re.replace_all(
        &content,
        format!(
            "// {{BEGIN_IMPL_PAR}}\n{}\n// {{END_IMPL_PAR}}",
            generated_code
        ),
    ).to_string();
    fs::write(lib_path, output).expect("Cannot write to lib.rs");
    println!("// Saved to lib.rs");
}
