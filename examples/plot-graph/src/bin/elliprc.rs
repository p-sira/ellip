/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::elliprc;
use ellip_plot_graph::*;
use plotly::{
    ImageFormat, Layout, Plot, Scatter,
    color::NamedColor,
    common::{Line, Mode},
    layout::{Annotation, Axis},
};

macro_rules! get_trace {
    ($x: expr, $y: expr, $name: expr) => {{
        let value = $x
            .iter()
            .map(|&xi| match elliprc(xi, $y) {
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
        get_trace!(&x, 1.0, "y=1", NamedColor::Red),
        get_trace!(&x, -1.0, "y=-1", NamedColor::Blue),
    ]);
    plot.set_layout(
        Layout::new()
            .title("Degenerate Symmetric Elliptic Integral of RF (RC)")
            .x_axis(Axis::new().title("x").show_line(true))
            .y_axis(
                Axis::new()
                    .title("RC(x,y)")
                    .show_line(true)
                    .range(vec![0.0, 2.0]),
            )
            .annotations(vec![Annotation::new()
            .text(format!(
                "Generated using <a href=\"https://docs.rs/ellip/latest/ellip/carlson/fn.elliprc.html\" target=\"_blank\">elliprc</a> from <a href=\"https://crates.io/crates/ellip\" target=\"_blank\">ellip</a> v{}",
                ellip_version()
            ))
                .x_ref("paper")
                .y_ref("paper")
                .y(-0.15)
                .x(1.08)
                .show_arrow(false)]),
    );

    make_html!(plot, "elliprc_plot.html");
    write_svg!(plot, "elliprc_plot.svg", 1000, 600, 1.0);
    println!("Done");
    Ok(())
}
