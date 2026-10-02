/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::{
    bulirsch::{cel, cel3},
    ellipe, ellipk, ellipke,
};
use ellip_dev_utils::{env::get_env, stats::Stats, test_report::err_func};
use std::time::Instant;

const PI: f64 = std::f64::consts::PI;

/// Fused Bulirsch complete elliptic evaluation for axial cylinder geometry (magba candidate).
#[inline]
fn cel_axial(kc: f64, gamma: f64) -> Result<(f64, f64), &'static str> {
    let mut kc = kc.abs();
    let mut aa_r: f64 = 0.0;
    let mut bb_r: f64 = -1.0;
    let mut c_r: f64 = 1.0;

    let pp_z: f64 = gamma.abs();
    let mut aa_z: f64 = 1.0;
    let mut bb_z: f64 = gamma / pp_z;
    let mut pp_z: f64 = pp_z;

    let mut e: f64 = kc;
    let mut m: f64 = 1.0;

    for _ in 0..16 {
        bb_r = (c_r * kc + bb_r) * 2.0;
        c_r = aa_r;

        let f_z = aa_z;
        let inv_pp_z = 1.0 / pp_z;
        aa_z = bb_z * inv_pp_z + aa_z;
        let g_z = e * inv_pp_z;
        bb_z = 2.0 * (f_z * g_z + bb_z);
        pp_z = g_z + pp_z;

        let m0 = m;
        m = kc + m;
        aa_r = bb_r / m + aa_r;

        if (m0 - kc).abs() > m0 * 1e-8 {
            kc = 2.0 * e.sqrt();
            e = kc * m;
            continue;
        }

        let ans_r = PI / 4.0 * aa_r / m;
        let ans_z = (PI / 2.0) * (aa_z * m + bb_z) / (m * (m + pp_z));
        return Ok((ans_r, ans_z));
    }

    Ok((
        cel(kc, 1.0, 1.0, -1.0)?,
        cel(kc, gamma * gamma, 1.0, gamma)?,
    ))
}

// =========================================================================
// Axial Cylinder B-field Kernels
// =========================================================================

#[inline]
fn unit_axial_cylinder_baseline(r: f64, z: f64, z0: f64) -> [f64; 3] {
    let (zp, zm) = (z + z0, z - z0);
    let (rp, rm) = (1.0 + r, 1.0 - r);

    let (zp2, zm2) = (zp * zp, zm * zm);
    let (rp2, rm2) = (rp * rp, rm * rm);

    let sq0 = f64::sqrt(zm2 + rp2);
    let sq1 = f64::sqrt(zp2 + rp2);

    let kp = f64::sqrt((zp2 + rm2) / (zp2 + rp2));
    let km = f64::sqrt((zm2 + rm2) / (zm2 + rp2));

    let gamma = rm / rp;
    let gamma2 = gamma * gamma;

    let br = (cel(kp, 1.0, 1.0, -1.0).unwrap() / sq1 - cel(km, 1.0, 1.0, -1.0).unwrap() / sq0) / PI;
    let bz = (zp * cel(kp, gamma2, 1.0, gamma).unwrap() / sq1
        - zm * cel(km, gamma2, 1.0, gamma).unwrap() / sq0)
        / (rp * PI);
    [br, 0.0, bz]
}

#[inline]
fn unit_axial_cylinder_specialized(r: f64, z: f64, z0: f64) -> [f64; 3] {
    let (zp, zm) = (z + z0, z - z0);
    let (rp, rm) = (1.0 + r, 1.0 - r);

    let (zp2, zm2) = (zp * zp, zm * zm);
    let (rp2, rm2) = (rp * rp, rm * rm);

    let sq0 = f64::sqrt(zm2 + rp2);
    let sq1 = f64::sqrt(zp2 + rp2);

    let kp = f64::sqrt((zp2 + rm2) / (zp2 + rp2));
    let km = f64::sqrt((zm2 + rm2) / (zm2 + rp2));

    let gamma = rm / rp;

    let (cr_p, cz_p) = cel_axial(kp, gamma).unwrap();
    let (cr_m, cz_m) = cel_axial(km, gamma).unwrap();

    let br = (cr_p / sq1 - cr_m / sq0) / PI;
    let bz = (zp * cz_p / sq1 - zm * cz_m / sq0) / (rp * PI);
    [br, 0.0, bz]
}

// =========================================================================
// Diametric Cylinder B-field Kernels
// =========================================================================

#[inline]
fn unit_diametric_cylinder_baseline(r: f64, phi: f64, z: f64, z0: f64) -> [f64; 3] {
    let (zp, zm) = (z + z0, z - z0);
    let (zp2, zm2) = (zp * zp, zm * zm);
    let r2 = r * r;

    let (rp, rm) = (r + 1.0, r - 1.0);
    let (_rp2, rm2) = (rp * rp, rm * rm);

    let (ap2, am2) = (zp2 + rm2, zm2 + rm2);
    let (ap, am) = (f64::sqrt(ap2), f64::sqrt(am2));

    let (argp, argm) = (-4.0 * r / ap2, -4.0 * r / am2);

    let (argc, one_over_rm) = if rm == 0.0 {
        (1e16, 0.0)
    } else {
        (-4.0 * r / rm2, 1.0 / rm)
    };

    let (ellk_p, ellk_m) = (ellipk(argp).unwrap(), ellipk(argm).unwrap());
    let (elle_p, elle_m) = (ellipe(argp).unwrap(), ellipe(argm).unwrap());
    let (ellpi_p, ellpi_m) = (
        cel(f64::sqrt(1.0 - argp), 1.0 - argc, 1.0, 1.0).unwrap(),
        cel(f64::sqrt(1.0 - argm), 1.0 - argc, 1.0, 1.0).unwrap(),
    );

    let br = -f64::cos(phi) / (4.0 * PI * r2)
        * (-zm * am * elle_m + zp * ap * elle_p + zm / am * (2.0 + zm2) * ellk_m
            - zp / ap * (2.0 + zp2) * ellk_p
            + zm * (rm2 - zm2) / am * ellpi_m
            - zp * (rm2 - zp2) / ap * ellpi_p);
    let bphi = -f64::sin(phi) / (4.0 * PI * r2)
        * (-zm * am * elle_m + zp * ap * elle_p + zm / am * (2.0 + zm2 + 2.0 * r2) * ellk_m
            - zp / ap * (2.0 + zp2 + 2.0 * r2) * ellk_p
            + zm * (rm2 - zm2) / am * ellpi_m
            - zp * (rm2 - zp2) / ap * ellpi_p);
    let bz = -f64::cos(phi) / (4.0 * PI * r)
        * (2.0 * am * elle_m - 2.0 * ap * elle_p + 2.0 / am * (1.0 + r2 + zm2) * ellk_m
            - 2.0 / ap * (1.0 + r2 + zp2) * ellk_p
            + 2.0 * rp * one_over_rm * (zm2 / am * ellpi_m - zp2 / ap * ellpi_p));

    [br, bphi, bz]
}

#[inline]
fn unit_diametric_cylinder_specialized(r: f64, phi: f64, z: f64, z0: f64) -> [f64; 3] {
    let (zp, zm) = (z + z0, z - z0);
    let (zp2, zm2) = (zp * zp, zm * zm);
    let r2 = r * r;

    let (rp, rm) = (r + 1.0, r - 1.0);
    let (_rp2, rm2) = (rp * rp, rm * rm);

    let (ap2, am2) = (zp2 + rm2, zm2 + rm2);
    let (ap, am) = (f64::sqrt(ap2), f64::sqrt(am2));

    let (argp, argm) = (-4.0 * r / ap2, -4.0 * r / am2);

    let (argc, one_over_rm) = if rm == 0.0 {
        (1e16, 0.0)
    } else {
        (-4.0 * r / rm2, 1.0 / rm)
    };

    let ((ellk_p, elle_p), (ellk_m, elle_m)) = (ellipke(argp).unwrap(), ellipke(argm).unwrap());
    let (ellpi_p, ellpi_m) = (
        cel3(f64::sqrt(1.0 - argp), 1.0 - argc).unwrap(),
        cel3(f64::sqrt(1.0 - argm), 1.0 - argc).unwrap(),
    );

    let br = -f64::cos(phi) / (4.0 * PI * r2)
        * (-zm * am * elle_m + zp * ap * elle_p + zm / am * (2.0 + zm2) * ellk_m
            - zp / ap * (2.0 + zp2) * ellk_p
            + zm * (rm2 - zm2) / am * ellpi_m
            - zp * (rm2 - zp2) / ap * ellpi_p);
    let bphi = -f64::sin(phi) / (4.0 * PI * r2)
        * (-zm * am * elle_m + zp * ap * elle_p + zm / am * (2.0 + zm2 + 2.0 * r2) * ellk_m
            - zp / ap * (2.0 + zp2 + 2.0 * r2) * ellk_p
            + zm * (rm2 - zm2) / am * ellpi_m
            - zp * (rm2 - zp2) / ap * ellpi_p);
    let bz = -f64::cos(phi) / (4.0 * PI * r)
        * (2.0 * am * elle_m - 2.0 * ap * elle_p + 2.0 / am * (1.0 + r2 + zm2) * ellk_m
            - 2.0 / ap * (1.0 + r2 + zp2) * ellk_p
            + 2.0 * rp * one_over_rm * (zm2 / am * ellpi_m - zp2 / ap * ellpi_p));

    [br, bphi, bz]
}

fn main() {
    let env = get_env();
    println!("===============================================================================");
    println!(" BENCHMARK REPORT: ellip 1.2.0 Specialized Cylinder Functions vs Baseline");
    println!(" CPU: {} ({:.2} GHz)", env.cpu, env.clock_speed);
    println!(" Platform: {} | rustc: {}", env.platform, env.rust_version);
    println!("===============================================================================\n");

    // -------------------------------------------------------------------------
    // 1. cel_axial vs separate cel calls
    // -------------------------------------------------------------------------
    println!("--- 1. cel_axial(kc, gamma) vs (cel(kc, 1.0, 1.0, -1.0), cel(kc, gamma^2, 1.0, gamma)) ---");
    let n_samples = 200_000;
    let mut kc_samples = Vec::with_capacity(n_samples);
    let mut gamma_samples = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let t = (i as f64 + 0.5) / n_samples as f64;
        kc_samples.push(1e-4 + (1.0 - 2e-4) * t);
        gamma_samples.push(-0.98 + 1.96 * t);
    }

    // Accuracy
    let mut cr_errs = Vec::with_capacity(n_samples);
    let mut cz_errs = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let kc = kc_samples[i];
        let gamma = gamma_samples[i];
        let (cr_spec, cz_spec) = cel_axial(kc, gamma).unwrap();
        let cr_base = cel(kc, 1.0, 1.0, -1.0).unwrap();
        let cz_base = cel(kc, gamma * gamma, 1.0, gamma).unwrap();
        cr_errs.push(err_func(cr_spec, cr_base));
        cz_errs.push(err_func(cz_spec, cz_base));
    }
    let stats_cr = Stats::from_vec(&cr_errs);
    let stats_cz = Stats::from_vec(&cz_errs);
    println!("Accuracy (units of eps):");
    println!(
        "  Cr: Mean = {:.2} eps, P99 = {:.2} eps, Max = {:.2} eps",
        stats_cr.mean, stats_cr.p99, stats_cr.max
    );
    println!(
        "  Cz: Mean = {:.2} eps, P99 = {:.2} eps, Max = {:.2} eps",
        stats_cz.mean, stats_cz.p99, stats_cz.max
    );

    // Timing
    let iters = 5;
    let mut base_time = f64::INFINITY;
    for _ in 0..iters {
        let start = Instant::now();
        let mut dummy = 0.0;
        for i in 0..n_samples {
            let kc = kc_samples[i];
            let gamma = gamma_samples[i];
            let cr = cel(kc, 1.0, 1.0, -1.0).unwrap();
            let cz = cel(kc, gamma * gamma, 1.0, gamma).unwrap();
            dummy += cr + cz;
        }
        std::hint::black_box(dummy);
        let elapsed = start.elapsed().as_secs_f64();
        if elapsed < base_time {
            base_time = elapsed;
        }
    }

    let mut spec_time = f64::INFINITY;
    for _ in 0..iters {
        let start = Instant::now();
        let mut dummy = 0.0;
        for i in 0..n_samples {
            let kc = kc_samples[i];
            let gamma = gamma_samples[i];
            let (cr, cz) = cel_axial(kc, gamma).unwrap();
            dummy += cr + cz;
        }
        std::hint::black_box(dummy);
        let elapsed = start.elapsed().as_secs_f64();
        if elapsed < spec_time {
            spec_time = elapsed;
        }
    }

    let base_ns = base_time * 1e9 / n_samples as f64;
    let spec_ns = spec_time * 1e9 / n_samples as f64;
    let speedup_cel_axial = base_ns / spec_ns;
    println!("Performance ({} samples):", n_samples);
    println!("  Baseline (2x cel): {:.2} ns/eval", base_ns);
    println!("  Specialized (cel_axial): {:.2} ns/eval", spec_ns);
    println!("  Speedup: {:.2}x\n", speedup_cel_axial);

    // -------------------------------------------------------------------------
    // 2. ellipke(m) vs (ellipk(m), ellipe(m)) on negative m
    // -------------------------------------------------------------------------
    println!("--- 2. ellipke(m) vs (ellipk(m), ellipe(m)) on negative m in [-1000, 0] ---");
    let mut neg_m_samples = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let t = (i as f64 + 0.5) / n_samples as f64;
        neg_m_samples.push(-1e-4 - 999.999 * t);
    }

    let mut k_errs = Vec::with_capacity(n_samples);
    let mut e_errs = Vec::with_capacity(n_samples);
    for &m in &neg_m_samples {
        let (k_spec, e_spec) = ellipke(m).unwrap();
        let k_base = ellipk(m).unwrap();
        let e_base = ellipe(m).unwrap();
        k_errs.push(err_func(k_spec, k_base));
        e_errs.push(err_func(e_spec, e_base));
    }
    let stats_k = Stats::from_vec(&k_errs);
    let stats_e = Stats::from_vec(&e_errs);
    println!("Accuracy vs ellipk / ellipe (units of eps):");
    println!(
        "  K: Mean = {:.2} eps, P99 = {:.2} eps, Max = {:.2} eps",
        stats_k.mean, stats_k.p99, stats_k.max
    );
    println!(
        "  E: Mean = {:.2} eps, P99 = {:.2} eps, Max = {:.2} eps",
        stats_e.mean, stats_e.p99, stats_e.max
    );

    let mut base_time = f64::INFINITY;
    for _ in 0..iters {
        let start = Instant::now();
        let mut dummy = 0.0;
        for &m in &neg_m_samples {
            let k = ellipk(m).unwrap();
            let e = ellipe(m).unwrap();
            dummy += k + e;
        }
        std::hint::black_box(dummy);
        let elapsed = start.elapsed().as_secs_f64();
        if elapsed < base_time {
            base_time = elapsed;
        }
    }

    let mut spec_time = f64::INFINITY;
    for _ in 0..iters {
        let start = Instant::now();
        let mut dummy = 0.0;
        for &m in &neg_m_samples {
            let (k, e) = ellipke(m).unwrap();
            dummy += k + e;
        }
        std::hint::black_box(dummy);
        let elapsed = start.elapsed().as_secs_f64();
        if elapsed < spec_time {
            spec_time = elapsed;
        }
    }

    let base_ns = base_time * 1e9 / n_samples as f64;
    let spec_ns = spec_time * 1e9 / n_samples as f64;
    let speedup_ellipke = base_ns / spec_ns;
    println!("Performance ({} samples):", n_samples);
    println!("  Baseline (ellipk + ellipe): {:.2} ns/eval", base_ns);
    println!("  Specialized (ellipke): {:.2} ns/eval", spec_ns);
    println!("  Speedup: {:.2}x\n", speedup_ellipke);

    // -------------------------------------------------------------------------
    // 3. cel3(kc, p) vs cel(kc, p, 1.0, 1.0)
    // -------------------------------------------------------------------------
    println!("--- 3. cel3(kc, p) vs cel(kc, p, 1.0, 1.0) ---");
    let mut p_samples = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let t = (i as f64 + 0.5) / n_samples as f64;
        p_samples.push(0.01 + 5.0 * t);
    }

    let mut cel3_errs = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let kc = kc_samples[i];
        let p = p_samples[i];
        let c_spec = cel3(kc, p).unwrap();
        let c_base = cel(kc, p, 1.0, 1.0).unwrap();
        cel3_errs.push(err_func(c_spec, c_base));
    }
    let stats_cel3 = Stats::from_vec(&cel3_errs);
    println!("Accuracy (units of eps):");
    println!(
        "  Mean = {:.2} eps, P99 = {:.2} eps, Max = {:.2} eps",
        stats_cel3.mean, stats_cel3.p99, stats_cel3.max
    );

    let mut base_time = f64::INFINITY;
    for _ in 0..iters {
        let start = Instant::now();
        let mut dummy = 0.0;
        for i in 0..n_samples {
            let kc = kc_samples[i];
            let p = p_samples[i];
            let c = cel(kc, p, 1.0, 1.0).unwrap();
            dummy += c;
        }
        std::hint::black_box(dummy);
        let elapsed = start.elapsed().as_secs_f64();
        if elapsed < base_time {
            base_time = elapsed;
        }
    }

    let mut spec_time = f64::INFINITY;
    for _ in 0..iters {
        let start = Instant::now();
        let mut dummy = 0.0;
        for i in 0..n_samples {
            let kc = kc_samples[i];
            let p = p_samples[i];
            let c = cel3(kc, p).unwrap();
            dummy += c;
        }
        std::hint::black_box(dummy);
        let elapsed = start.elapsed().as_secs_f64();
        if elapsed < spec_time {
            spec_time = elapsed;
        }
    }

    let base_ns = base_time * 1e9 / n_samples as f64;
    let spec_ns = spec_time * 1e9 / n_samples as f64;
    let speedup_cel3 = base_ns / spec_ns;
    println!("Performance ({} samples):", n_samples);
    println!("  Baseline (cel): {:.2} ns/eval", base_ns);
    println!("  Specialized (cel3): {:.2} ns/eval", spec_ns);
    println!("  Speedup: {:.2}x\n", speedup_cel3);

    // -------------------------------------------------------------------------
    // 4. End-to-End Axial Cylinder Field Benchmark
    // -------------------------------------------------------------------------
    println!("--- 4. End-to-End unit_axial_cylinder_B_cyl ---");
    let grid_r = 300;
    let grid_z = 300;
    let total_points = grid_r * grid_z;
    let mut grid_points = Vec::with_capacity(total_points);
    for ir in 0..grid_r {
        let r = 0.06 + 4.0 * (ir as f64 / (grid_r - 1) as f64);
        for iz in 0..grid_z {
            let z = -3.0 + 6.0 * (iz as f64 / (grid_z - 1) as f64);
            grid_points.push((r, z));
        }
    }
    let z0 = 1.0;

    let mut br_errs = Vec::with_capacity(total_points);
    let mut bz_errs = Vec::with_capacity(total_points);
    for &(r, z) in &grid_points {
        let base = unit_axial_cylinder_baseline(r, z, z0);
        let spec = unit_axial_cylinder_specialized(r, z, z0);
        br_errs.push(err_func(spec[0], base[0]));
        bz_errs.push(err_func(spec[2], base[2]));
    }
    let stats_ax_br = Stats::from_vec(&br_errs);
    let stats_ax_bz = Stats::from_vec(&bz_errs);
    println!(
        "Accuracy across {} grid points (units of eps):",
        total_points
    );
    println!(
        "  Br: Mean = {:.2} eps, P99 = {:.2} eps, Max = {:.2} eps",
        stats_ax_br.mean, stats_ax_br.p99, stats_ax_br.max
    );
    println!(
        "  Bz: Mean = {:.2} eps, P99 = {:.2} eps, Max = {:.2} eps",
        stats_ax_bz.mean, stats_ax_bz.p99, stats_ax_bz.max
    );

    let mut base_time = f64::INFINITY;
    for _ in 0..iters {
        let start = Instant::now();
        let mut dummy = 0.0;
        for &(r, z) in &grid_points {
            let b = unit_axial_cylinder_baseline(r, z, z0);
            dummy += b[0] + b[2];
        }
        std::hint::black_box(dummy);
        let elapsed = start.elapsed().as_secs_f64();
        if elapsed < base_time {
            base_time = elapsed;
        }
    }

    let mut spec_time = f64::INFINITY;
    for _ in 0..iters {
        let start = Instant::now();
        let mut dummy = 0.0;
        for &(r, z) in &grid_points {
            let b = unit_axial_cylinder_specialized(r, z, z0);
            dummy += b[0] + b[2];
        }
        std::hint::black_box(dummy);
        let elapsed = start.elapsed().as_secs_f64();
        if elapsed < spec_time {
            spec_time = elapsed;
        }
    }

    let base_ns = base_time * 1e9 / total_points as f64;
    let spec_ns = spec_time * 1e9 / total_points as f64;
    let speedup_axial = base_ns / spec_ns;
    println!("Performance ({} field evaluations):", total_points);
    println!(
        "  Baseline: {:.2} ns/eval ({:.2} ms total)",
        base_ns,
        base_time * 1e3
    );
    println!(
        "  Specialized: {:.2} ns/eval ({:.2} ms total)",
        spec_ns,
        spec_time * 1e3
    );
    println!("  Speedup: {:.2}x\n", speedup_axial);

    // -------------------------------------------------------------------------
    // 5. End-to-End Diametric Cylinder Field Benchmark
    // -------------------------------------------------------------------------
    println!("--- 5. End-to-End unit_diametric_cylinder_B_cyl ---");
    let phi = 0.7853981633974483; // pi / 4
    let mut br_dia_errs = Vec::with_capacity(total_points);
    let mut bphi_dia_errs = Vec::with_capacity(total_points);
    let mut bz_dia_errs = Vec::with_capacity(total_points);
    for &(r, z) in &grid_points {
        let base = unit_diametric_cylinder_baseline(r, phi, z, z0);
        let spec = unit_diametric_cylinder_specialized(r, phi, z, z0);
        br_dia_errs.push(err_func(spec[0], base[0]));
        bphi_dia_errs.push(err_func(spec[1], base[1]));
        bz_dia_errs.push(err_func(spec[2], base[2]));
    }
    let stats_dia_br = Stats::from_vec(&br_dia_errs);
    let stats_dia_bphi = Stats::from_vec(&bphi_dia_errs);
    let stats_dia_bz = Stats::from_vec(&bz_dia_errs);
    println!(
        "Accuracy across {} grid points (units of eps):",
        total_points
    );
    println!(
        "  Br:   Mean = {:.2} eps, P99 = {:.2} eps, Max = {:.2} eps",
        stats_dia_br.mean, stats_dia_br.p99, stats_dia_br.max
    );
    println!(
        "  Bphi: Mean = {:.2} eps, P99 = {:.2} eps, Max = {:.2} eps",
        stats_dia_bphi.mean, stats_dia_bphi.p99, stats_dia_bphi.max
    );
    println!(
        "  Bz:   Mean = {:.2} eps, P99 = {:.2} eps, Max = {:.2} eps",
        stats_dia_bz.mean, stats_dia_bz.p99, stats_dia_bz.max
    );

    let mut base_time = f64::INFINITY;
    for _ in 0..iters {
        let start = Instant::now();
        let mut dummy = 0.0;
        for &(r, z) in &grid_points {
            let b = unit_diametric_cylinder_baseline(r, phi, z, z0);
            dummy += b[0] + b[1] + b[2];
        }
        std::hint::black_box(dummy);
        let elapsed = start.elapsed().as_secs_f64();
        if elapsed < base_time {
            base_time = elapsed;
        }
    }

    let mut spec_time = f64::INFINITY;
    for _ in 0..iters {
        let start = Instant::now();
        let mut dummy = 0.0;
        for &(r, z) in &grid_points {
            let b = unit_diametric_cylinder_specialized(r, phi, z, z0);
            dummy += b[0] + b[1] + b[2];
        }
        std::hint::black_box(dummy);
        let elapsed = start.elapsed().as_secs_f64();
        if elapsed < spec_time {
            spec_time = elapsed;
        }
    }

    let base_ns = base_time * 1e9 / total_points as f64;
    let spec_ns = spec_time * 1e9 / total_points as f64;
    let speedup_diametric = base_ns / spec_ns;
    println!("Performance ({} field evaluations):", total_points);
    println!(
        "  Baseline: {:.2} ns/eval ({:.2} ms total)",
        base_ns,
        base_time * 1e3
    );
    println!(
        "  Specialized: {:.2} ns/eval ({:.2} ms total)",
        spec_ns,
        spec_time * 1e3
    );
    println!("  Speedup: {:.2}x\n", speedup_diametric);

    println!("===============================================================================");
    println!(" SUMMARY OF SPEEDUPS:");
    println!(
        "  cel_axial vs 2x cel:               {:.2}x",
        speedup_cel_axial
    );
    println!(
        "  ellipke vs (ellipk + ellipe):      {:.2}x",
        speedup_ellipke
    );
    println!("  cel3 vs cel:                       {:.2}x", speedup_cel3);
    println!("  End-to-end Axial Cylinder B:       {:.2}x", speedup_axial);
    println!(
        "  End-to-end Diametric Cylinder B:   {:.2}x",
        speedup_diametric
    );
    println!("===============================================================================");
}
