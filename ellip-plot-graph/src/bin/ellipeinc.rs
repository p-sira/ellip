/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use std::f64::consts::{FRAC_PI_2, FRAC_PI_3, FRAC_PI_4};

use ellip::ellipeinc;
use ellip_plot_graph::*;
use plotly::{
    Layout, Plot, Scatter,
    color::NamedColor,
    common::{Line, Mode},
    layout::{Annotation, Axis},
};

macro_rules! get_trace {
    ($phi: expr, $m: expr, $name: expr) => {{
        let value = $m
            .iter()
            .map(|&mi| match ellipeinc($phi, mi) {
                Ok(ans) => ans,
                Err(_) => f64::NAN,
            })
            .collect();

        Scatter::new($m.clone(), value)
            .mode(Mode::Lines)
            .name($name)
    }};
    ($phi: expr, $m: expr, $name: expr, $line_color: expr) => {
        get_trace!($phi, $m, $name).line(Line::new().color($line_color))
    };
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n_points = 100;
    let range_m = [-5, 2];
    let m: Vec<f64> = (range_m[0] * n_points..=range_m[1] * n_points)
        .map(|x| x as f64 / n_points as f64)
        .collect();

    let mut plot = Plot::new();
    plot.add_traces(vec![
        get_trace!(FRAC_PI_2, &m, "φ=π/2", NamedColor::Red),
        get_trace!(FRAC_PI_3, &m, "φ=π/3"),
        get_trace!(FRAC_PI_4, &m, "φ=π/4"),
    ]);
    plot.set_layout(
        Layout::new()
            .title("Incomplete Elliptic Integral of the First Kind (E)")
            .x_axis(Axis::new().title("m").show_line(true))
            .y_axis(
                Axis::new()
                    .title("E(φ,m)")
                    .show_line(true)
                    .range(vec![0.0, 4.0]),
            )
            .legend(
                plotly::layout::Legend::new()
                    .x(1.0)
                    .x_anchor(plotly::common::Anchor::Right)
            )
            .annotations(vec![Annotation::new()
                .text(format!(
                    "Generated using <a href=\"https://docs.rs/ellip/latest/ellip/legendre/fn.ellipeinc.html\" target=\"_blank\">ellipeinc</a> from <a href=\"https://crates.io/crates/ellip\" target=\"_blank\">ellip</a> v{}",
                    ellip_version()
                ))
                .x_ref("paper")
                .y_ref("paper")
                .y(-0.15)
                .x(1.08)
                .show_arrow(false)]),
    );

    make_html!(plot, "ellipeinc.html");
    write_svg!(plot, "ellipeinc.svg", 1000, 600, 1.0);
    Ok(())
}
