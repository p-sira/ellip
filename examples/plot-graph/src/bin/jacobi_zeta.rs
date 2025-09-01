/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use std::f64::consts::PI;

use ellip::jacobi_zeta;
use ellip_plot_graph::*;
use plotly::{
    color::NamedColor,
    common::{Line, Mode},
    layout::{Annotation, Axis},
    ImageFormat, Layout, Plot, Scatter,
};

macro_rules! get_trace {
    ($phi: expr, $m: expr, $name: expr) => {{
        let value = $phi
            .iter()
            .map(|&phi_i| match jacobi_zeta(phi_i, $m) {
                Ok(ans) => ans,
                Err(_) => f64::NAN,
            })
            .collect();

        Scatter::new($phi.clone(), value)
            .mode(Mode::Lines)
            .name($name)
    }};
    ($phi: expr, $m: expr, $name: expr, $line_color: expr) => {
        get_trace!($phi, $m, $name).line(Line::new().color($line_color))
    };
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n_points = 1000;
    let range_phi = [-2.0 * PI, 2.0 * PI];
    let phi: Vec<f64> = (0..=n_points)
        .map(|i| range_phi[0] + i as f64 * (range_phi[1] - range_phi[0]) / n_points as f64)
        .collect();

    // Create tick values in units of π
    let tick_vals: Vec<f64> = (-2..=2).map(|i| i as f64 * PI).collect();
    let tick_text: Vec<String> = (-2..=2)
        .map(|i| {
            if i == 0 {
                "0".to_string()
            } else if i == 1 {
                "π".to_string()
            } else {
                format!("{}π", i)
            }
        })
        .collect();

    let mut plot = Plot::new();
    plot.add_traces(vec![
        get_trace!(&phi, 0.4, "m=0.4", NamedColor::Red),
        get_trace!(&phi, 0.7, "m=0.7", NamedColor::Blue),
        get_trace!(&phi, 0.99, "m=0.99", NamedColor::Green),
        get_trace!(&phi, 0.999999, "m=0.999999", NamedColor::Purple),
    ]);
    plot.set_layout(
        Layout::new()
            .title("Jacobi Zeta Function Z(φ,m)")
            .x_axis(
                Axis::new()
                    .title("φ")
                    .show_line(true)
                    .tick_values(tick_vals)
                    .tick_text(tick_text)
            )
            .y_axis(
                Axis::new()
                    .title("Z(φ,m)")
                    .show_line(true)
                    .range(vec![-1.0, 1.0]),
            )
            .annotations(vec![Annotation::new()
            .text(format!(
                "Generated using <a href=\"https://docs.rs/ellip/latest/ellip/fn.jacobi_zeta.html\" target=\"_blank\">jacobi_zeta</a> from <a href=\"https://crates.io/crates/ellip\" target=\"_blank\">ellip</a> v{}",
                ellip_version()
            ))
                .x_ref("paper")
                .y_ref("paper")
                .y(-0.15)
                .x(1.08)
                .show_arrow(false)]),
    );

    make_html!(plot, "jacobi_zeta_plot.html");
    write_svg!(plot, "jacobi_zeta_plot.svg", 900, 600, 1.0);
    println!("Done");
    Ok(())
}
