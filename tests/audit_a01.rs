#![cfg(not(feature = "test_force_fail"))]
// Regression for audit finding A1; references use exact binary inputs at 100+ digits.
use ellip::*;
#[allow(dead_code)]
fn close(actual: f64, expected: f64, rtol: f64) {
    assert!(actual.is_finite(), "actual={actual}, expected={expected}");
    assert!(
        (actual - expected).abs() <= rtol * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}
#[test]
fn amplitude_periods() {
    close(
        ellippiinc_bulirsch(2.0, 0.5, 0.5).unwrap(),
        3.8198568874384073,
        2e-15,
    );
    for phi in [std::f64::consts::PI, -std::f64::consts::PI, 7.0, -7.0] {
        close(ellippiinc_bulirsch(phi, 0.0, 0.0).unwrap(), phi, 2e-15);
        close(
            ellippiinc_bulirsch(phi, 0.5, 0.5).unwrap(),
            ellippiinc(phi, 0.5, 0.5).unwrap(),
            2e-15,
        );
    }
    for phi in [std::f32::consts::FRAC_PI_2, -std::f32::consts::FRAC_PI_2] {
        close(
            ellippiinc_bulirsch(phi, 0.5, 0.5).unwrap() as f64,
            (phi.signum() as f64) * 2.7012878857298321,
            3e-7,
        );
    }
}
