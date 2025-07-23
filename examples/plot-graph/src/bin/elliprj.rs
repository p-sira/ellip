/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::elliprj;
use ellip_plot_graph::*;
use plotly::{
    ImageFormat, Layout, Plot, Scatter,
    color::NamedColor,
    common::{Anchor, Line, Mode},
    layout::{Annotation, Axis, Legend},
};

macro_rules! get_trace {
    ($x: expr, $y: expr, $p: expr, $name: expr) => {{
        let value = $x
            .iter()
            .map(|&xi| match elliprj(xi, $y, 1.0, $p) {
                Ok(ans) => ans,
                Err(_) => f64::NAN,
            })
            .collect();

        Scatter::new($x.clone(), value)
            .mode(Mode::Lines)
            .name($name)
    }};
    ($x: expr, $y: expr, $p: expr, $name: expr, $line_color: expr) => {
        get_trace!($x, $y, $p, $name).line(Line::new().color($line_color))
    };
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n_points = 100;

    // Make the plot more dense near x=0
    let x_near_zero: Vec<f64> = (0..n_points)
        .map(|i| i as f64 / n_points as f64 / 10.0)
        .collect();

    let x_full: Vec<f64> = (0..n_points)
        .map(|i| 0.1 + 1e-5 + (1.0 - (0.1 + 1e-5)) * i as f64 / n_points as f64)
        .collect();

    let x = [x_near_zero, x_full].concat();

    let mut plot = Plot::new();
    plot.add_traces(vec![
        get_trace!(&x, 0.0, 2.0, "y=0.0, p=2", NamedColor::Red)
            .legend_group("p=2")
            .legend_group_title("p=2"),
        get_trace!(&x, 0.1, 2.0, "y=0.1, p=2").legend_group("p=2"),
        get_trace!(&x, 0.5, 2.0, "y=0.5, p=2").legend_group("p=2"),
        get_trace!(&x, 1.0, 2.0, "y=1.0, p=2", NamedColor::Blue).legend_group("p=2"),
        get_trace!(&x, 1.0, -2.0, "y=1.0, p=-2", NamedColor::Navy)
            .legend_group("p=-2")
            .legend_group_title("p=-2"),
        get_trace!(&x, 0.5, -2.0, "y=0.5, p=-2", NamedColor::ForestGreen).legend_group("p=-2"),
        get_trace!(&x, 0.1, -2.0, "y=0.1, p=-2", NamedColor::DarkOrange).legend_group("p=-2"),
        get_trace!(&x, 0.0, -2.0, "y=0.0, p=-2", NamedColor::Crimson).legend_group("p=-2"),
    ]);
    plot.set_layout(
        Layout::new()
            .title("Symmetric Elliptic Integral of the Third Kind (RJ)")
            .x_axis(Axis::new().title("x").show_line(true))
            .y_axis(
                Axis::new()
                    .title("RJ(x,y,1,p)")
                    .show_line(true)
                    .range(vec![-5.0, 5.0]),
            ).legend(Legend::new().y_anchor(Anchor::Middle).y(0.5))
            .annotations(vec![Annotation::new()
            .text(format!(
                "Generated using the function <a href=\"https://docs.rs/ellip/latest/ellip/carlson/fn.elliprj.html\" target=\"_blank\">elliprj</a> from <a href=\"https://crates.io/crates/ellip\" target=\"_blank\">ellip</a> v{}",
                env!("CARGO_PKG_VERSION")
            ))
                .x_ref("paper")
                .y_ref("paper")
                .y(-0.15)
                .x(1.08)
                .show_arrow(false)]),
    );

    make_html!("elliprj_plot.html");
    plot.write_image(
        figure_path!("elliprj_plot.svg"),
        ImageFormat::SVG,
        900,
        600,
        1.0,
    );
    Ok(())
}
