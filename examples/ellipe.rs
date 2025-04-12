/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::ellipe;
use plotly::{
    color::NamedColor,
    common::{Line, Mode},
    layout::{Annotation, Axis},
    ImageFormat, Layout, Plot, Scatter,
};


fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n_points = 100;
    let range_m = [-2, 1];

    let m: Vec<f64> = (range_m[0] * n_points..range_m[1] * n_points)
        .map(|x| x as f64 / n_points as f64)
        .collect();

    let ellipe_values: Vec<f64> = m.iter().map(|&mi| ellipe(mi).unwrap()).collect();

    let trace = Scatter::new(m, ellipe_values)
        .mode(Mode::Lines)
        .name("E(m)")
        .line(Line::new().color(NamedColor::Red));

    let mut plot = Plot::new();
    plot.add_trace(trace);
    plot.set_layout(
        Layout::new()
            .title("Complete Elliptic Integral of the Second Kind (E)")
            .x_axis(Axis::new().title("m").show_line(true))
            .y_axis(
                Axis::new()
                    .title("E(m)")
                    .show_line(true)
                    .range(vec![0.0, 4.0]),
            )
            .annotations(vec![Annotation::new()
            .text(format!(
                "Generated using the function <a href=\"https://docs.rs/ellip/latest/ellip/legendre/fn.ellipe.html\" target=\"_blank\">ellipe</a> from <a href=\"https://crates.io/crates/ellip\" target=\"_blank\">ellip</a> v{}",
                env!("CARGO_PKG_VERSION")
            ))
            .x_ref("paper")
            .y_ref("paper")
            .y(-0.15)
            .x(1.05)
            .show_arrow(false)]),
    );

    plot.show_html("figures/ellipe_plot.html");
    plot.write_image("figures/ellipe_plot.svg", ImageFormat::SVG, 900, 600, 1.0);
    Ok(())
}
