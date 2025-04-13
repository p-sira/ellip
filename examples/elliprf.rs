/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::elliprf;
use plotly::{
    color::NamedColor,
    common::{Line, Mode},
    layout::{Annotation, Axis},
    ImageFormat, Layout, Plot, Scatter,
};

macro_rules! get_trace {
    ($x: expr, $y: expr, $name: expr) => {{
        let value = $x
            .iter()
            .map(|&xi| match elliprf(xi, $y, 1.0) {
                Ok(ans) => ans,
                Err(_) => f64::NAN,
            })
            .collect();

        Scatter::new($x.clone(), value)
            .mode(Mode::Lines)
            .name($name)
    }};
    ($x: expr, $y: expr, $name: expr, $line_color: expr) => {
        get_trace!($x, $y, $name).line(Line::new().color($line_color))
    };
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n_points = 100;
    let range_x = [0, 5];
    let x: Vec<f64> = (range_x[0] * n_points..=range_x[1] * n_points)
        .map(|x| x as f64 / n_points as f64)
        .collect();

    let mut plot = Plot::new();
    plot.add_traces(vec![
        get_trace!(&x, 0.0, "y=0", NamedColor::Red),
        get_trace!(&x, 0.1, "y=0.1"),
        get_trace!(&x, 0.5, "y=0.5"),
        get_trace!(&x, 1.0, "y=1.0", NamedColor::Blue),
    ]);
    plot.set_layout(
        Layout::new()
            .title("Symmetric Elliptic Integral of the First Kind (RF)")
            .x_axis(Axis::new().title("x").show_line(true))
            .y_axis(
                Axis::new()
                    .title("RF(x,y,1)")
                    .show_line(true)
                    .range(vec![0.0, 3.0]),
            )
            .annotations(vec![Annotation::new()
            .text(format!(
                "Generated using the function <a href=\"https://docs.rs/ellip/latest/ellip/carlson/fn.elliprf.html\" target=\"_blank\">elliprf</a> from <a href=\"https://crates.io/crates/ellip\" target=\"_blank\">ellip</a> v{}",
                env!("CARGO_PKG_VERSION")
            ))
                .x_ref("paper")
                .y_ref("paper")
                .y(-0.15)
                .x(1.08)
                .show_arrow(false)]),
    );

    plot.show_html("figures/elliprf_plot.html");
    plot.write_image("figures/elliprf_plot.svg", ImageFormat::SVG, 900, 600, 1.0);
    Ok(())
}
