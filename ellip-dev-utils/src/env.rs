/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

pub struct Env {
    pub rust_version: String,
    pub platform: String,
    pub ellip_version: String,
    pub cpu: String,
    /// Clock speed in Hz
    pub clock_speed: u64,
    pub total_memory: u64,
}

pub fn get_env() -> Env {
    use std::fs;
    use std::process::Command;

    let clock_speed = {
        // sysfs cpuinfo_max_freq (kHz)
        fs::read_to_string("/sys/devices/system/cpu/cpu0/cpufreq/cpuinfo_max_freq")
            .ok()
            .and_then(|s| s.trim().parse::<u64>().ok())
            .map(|khz| khz * 1000) //Hz
            .unwrap()
    };

    let mut sys = sysinfo::System::new_all();
    sys.refresh_all();
    let cpu = match sys.cpus().first() {
        Some(cpu) => cpu.brand(),
        None => "Unknown CPU",
    }
    .to_string();

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

    let total_memory = sys.total_memory();

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

    Env {
        rust_version,
        platform,
        ellip_version,
        cpu,
        clock_speed,
        total_memory,
    }
}

pub fn format_clock_speed(clock_speed: u64) -> String {
    if clock_speed > 10_000_000 {
        format!("{} GHz", clock_speed / 1_000_000_000)
    } else {
        format!("{} MHz", clock_speed / 1_000_000)
    }
}
