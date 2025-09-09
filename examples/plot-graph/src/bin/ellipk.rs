/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::ellipk;
use ellip_plot_graph::*;
use plotly::{
    ImageFormat, Layout, Plot, Scatter,
    color::NamedColor,
    common::{Line, Mode},
    layout::{Annotation, Axis},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n_points = 100;
    let range_m = [-2, 1];

    let m: Vec<f64> = (range_m[0] * n_points..=range_m[1] * n_points)
        .map(|x| x as f64 / n_points as f64)
        .collect();

    let ellipk_values: Vec<f64> = m.iter().map(|&mi| ellipk(mi).unwrap()).collect();

    let trace = Scatter::new(m, ellipk_values)
        .mode(Mode::Lines)
        .name("K(m)")
        .line(Line::new().color(NamedColor::Red));

    let mut plot = Plot::new();
    plot.add_trace(trace);
    plot.set_layout(
        Layout::new()
            .title("Complete Elliptic Integral of the First Kind (K)")
            .x_axis(Axis::new().title("m").show_line(true))
            .y_axis(
                Axis::new()
                    .title("K(m)")
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
                "Generated using <a href=\"https://docs.rs/ellip/latest/ellip/legendre/fn.ellipk.html\" target=\"_blank\">ellipk</a> from <a href=\"https://crates.io/crates/ellip\" target=\"_blank\">ellip</a> v{}",
                ellip_version()
            ))
            .x_ref("paper")
            .y_ref("paper")
            .y(-0.15)
            .x(1.05)
            .show_arrow(false)]),
    );

    make_html!(plot, "ellipk_plot.html");
    write_svg!(plot, "ellipk_plot.svg", 1000, 600, 1.0);
    println!("Done");
    Ok(())
}
