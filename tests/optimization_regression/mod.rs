use ellip::{ellipd, ellipdinc, ellipe, ellipeinc, ellipf, ellipk};

#[test]
fn principal_angle_and_period_boundaries_f64() {
    let half_pi = core::f64::consts::FRAC_PI_2;
    let phi_below = f64::from_bits(half_pi.to_bits() - 1);
    let phi_above = f64::from_bits(half_pi.to_bits() + 1);
    let angles = [
        f64::MIN_POSITIVE,
        1e-12,
        0.2,
        0.9,
        phi_below,
        half_pi,
        phi_above,
        core::f64::consts::PI,
        3.0 * half_pi,
        9.0 * half_pi,
        2.0_f64.powi(40) * core::f64::consts::PI,
        2.0_f64.powi(49) * core::f64::consts::PI,
    ];
    for m in [-2.0, -0.4, 0.5, 0.99] {
        let complete = (ellipk(m).unwrap(), ellipe(m).unwrap(), ellipd(m).unwrap());
        for phi in angles {
            let actual = (
                ellipf(phi, m).unwrap(),
                ellipeinc(phi, m).unwrap(),
                ellipdinc(phi, m).unwrap(),
            );
            let negated = (
                ellipf(-phi, m).unwrap(),
                ellipeinc(-phi, m).unwrap(),
                ellipdinc(-phi, m).unwrap(),
            );
            let next_period = (
                ellipf(phi + core::f64::consts::PI, m).unwrap(),
                ellipeinc(phi + core::f64::consts::PI, m).unwrap(),
                ellipdinc(phi + core::f64::consts::PI, m).unwrap(),
            );
            for (value, negative, shifted, period) in [
                (actual.0, negated.0, next_period.0, complete.0),
                (actual.1, negated.1, next_period.1, complete.1),
                (actual.2, negated.2, next_period.2, complete.2),
            ] {
                assert!(
                    (value + negative).abs() <= 1e-14 * value.abs().max(1.0),
                    "oddness: phi={phi}, m={m}, value={value}, negative={negative}"
                );
                let expected = value + 2.0 * period;
                assert!(
                    (shifted - expected).abs() <= 1e-14 * expected.abs().max(1.0),
                    "period: phi={phi}, m={m}, shifted={shifted}, expected={expected}"
                );
            }
        }
    }
}

#[test]
fn principal_angle_and_period_boundaries_f32() {
    let half_pi = core::f32::consts::FRAC_PI_2;
    let angles = [
        f32::MIN_POSITIVE,
        0.2,
        f32::from_bits(half_pi.to_bits() - 1),
        half_pi,
        f32::from_bits(half_pi.to_bits() + 1),
        3.0 * half_pi,
        2.0_f32.powi(20) * core::f32::consts::PI,
    ];
    for m in [-2.0_f32, -0.4, 0.5, 0.99] {
        let complete = (ellipk(m).unwrap(), ellipe(m).unwrap(), ellipd(m).unwrap());
        for phi in angles {
            let actual = (
                ellipf(phi, m).unwrap(),
                ellipeinc(phi, m).unwrap(),
                ellipdinc(phi, m).unwrap(),
            );
            let negative = (
                ellipf(-phi, m).unwrap(),
                ellipeinc(-phi, m).unwrap(),
                ellipdinc(-phi, m).unwrap(),
            );
            let next = (
                ellipf(phi + core::f32::consts::PI, m).unwrap(),
                ellipeinc(phi + core::f32::consts::PI, m).unwrap(),
                ellipdinc(phi + core::f32::consts::PI, m).unwrap(),
            );
            for (value, negated, shifted, period) in [
                (actual.0, negative.0, next.0, complete.0),
                (actual.1, negative.1, next.1, complete.1),
                (actual.2, negative.2, next.2, complete.2),
            ] {
                assert!(
                    (value + negated).abs() <= 3e-6 * value.abs().max(1.0),
                    "oddness: phi={phi}, m={m}"
                );
                let expected = value + 2.0 * period;
                assert!(
                    (shifted - expected).abs() <= 3e-6 * expected.abs().max(1.0),
                    "period: phi={phi}, m={m}, shifted={shifted}, expected={expected}"
                );
            }
        }
    }
}

#[test]
fn invalid_amplitudes_do_not_panic_in_period_reduction() {
    for m in [0.25_f64, -2.0] {
        assert!(ellipf(f64::NAN, m).is_err());
        assert!(ellipeinc(f64::NAN, m).is_err());
        assert!(ellipdinc(f64::NAN, m).is_err());
    }
}
