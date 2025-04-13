/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::cel1;
use plotly::{
    color::NamedColor,
    common::{Line, Mode},
    layout::{Annotation, Axis},
    ImageFormat, Layout, Plot, Scatter,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n_points = 100;
    let range_kc = [-3, 3];

    let kc: Vec<f64> = (range_kc[0] * n_points..=range_kc[1] * n_points)
        .map(|x| x as f64 / n_points as f64)
        .collect();

    let cel1_values: Vec<f64> = kc
        .iter()
        .map(|&kci| match cel1(kci) {
            Ok(ans) => ans,
            Err(_) => f64::NAN,
        })
        .collect();

    let trace = Scatter::new(kc, cel1_values)
        .mode(Mode::Lines)
        .name("cel1(kc)")
        .line(Line::new().color(NamedColor::Red));

    let mut plot = Plot::new();
    plot.add_trace(trace);
    plot.set_layout(
        Layout::new()
            .title("Bulirsch's Complete Elliptic Integral of the First Kind (cel1)")
            .x_axis(Axis::new().title("kc").show_line(true))
            .y_axis(
                Axis::new()
                    .title("cel1(kc)")
                    .show_line(true)
                    .range(vec![0.0, 6.0]),
            )
            .annotations(vec![Annotation::new()
            .text(format!(
                "Generated using the function <a href=\"https://docs.rs/ellip/latest/ellip/bulirsch/fn.cel1.html\" target=\"_blank\">cel1</a> from <a href=\"https://crates.io/crates/ellip\" target=\"_blank\">ellip</a> v{}",
                env!("CARGO_PKG_VERSION")
            ))
            .x_ref("paper")
            .y_ref("paper")
            .y(-0.15)
            .x(1.05)
            .show_arrow(false)]),
    );

    plot.show_html("figures/cel1_plot.html");
    plot.write_image("figures/cel1_plot.svg", ImageFormat::SVG, 900, 600, 1.0);
    Ok(())
}
