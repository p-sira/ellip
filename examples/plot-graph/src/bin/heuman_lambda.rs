/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use std::f64::consts::PI;

use ellip::heuman_lambda;
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
            .map(|&phi_i| match heuman_lambda(phi_i, $m) {
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
    let range_phi = [-PI, PI];
    let phi: Vec<f64> = (0..=n_points)
        .map(|i| range_phi[0] + i as f64 * (range_phi[1] - range_phi[0]) / n_points as f64)
        .collect();

    let mut plot = Plot::new();
    plot.add_traces(vec![
        get_trace!(&phi, 0.1, "m=0.1", NamedColor::Red),
        get_trace!(&phi, 0.8, "m=0.8", NamedColor::Blue),
        get_trace!(&phi, 0.99, "m=0.99", NamedColor::Green),
    ]);
    plot.set_layout(
        Layout::new()
            .title("Heuman Lambda Function Λ0(φ,m)")
            .x_axis(
                Axis::new()
                    .title("φ")
                    .show_line(true)
            )
            .y_axis(
                Axis::new()
                    .title("Λ0(φ,m)")
                    .show_line(true)
                    .range(vec![-3.0, 3.0]),
            )
            .annotations(vec![Annotation::new()
            .text(format!(
                "Generated using <a href=\"https://docs.rs/ellip/latest/ellip/fn.heuman_lambda.html\" target=\"_blank\">heuman_lambda</a> from <a href=\"https://crates.io/crates/ellip\" target=\"_blank\">ellip</a> v{}",
                ellip_version()
            ))
                .x_ref("paper")
                .y_ref("paper")
                .y(-0.15)
                .x(1.08)
                .show_arrow(false)]),
    );

    make_html!(plot, "heuman_lambda_plot.html");
    write_svg!(plot, "heuman_lambda_plot.svg", 900, 600, 1.0);
    println!("Done");
    Ok(())
}
