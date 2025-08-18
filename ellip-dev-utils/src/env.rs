/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

pub fn get_env() -> [String; 3] {
    use std::process::Command;

    let rust_version = {
        let output = Command::new("rustc")
            .arg("--version")
            .output()
            .unwrap()
            .stdout;
        String::from_utf8_lossy(&output)
            .split_whitespace()
            .collect::<Vec<&str>>()[1]
            .to_owned()
    };

    let platform = {
        let output = Command::new("rustc").arg("-vV").output().unwrap().stdout;
        String::from_utf8_lossy(&output)
            .lines()
            .find(|line| line.starts_with("host:"))
            .unwrap()
            .strip_prefix("host: ")
            .unwrap()
            .to_owned()
    };

    let ellip_version: String = {
        let output = Command::new("cargo")
            .args(["tree", "--invert", "--package", "ellip"])
            .output()
            .unwrap()
            .stdout;
        String::from_utf8_lossy(&output)
            .lines()
            .next()
            .and_then(|line| line.strip_prefix("ellip v"))
            .unwrap()
            .split_whitespace()
            .next()
            .unwrap()
            .to_owned()
    };

    [rust_version, platform, ellip_version]
}
