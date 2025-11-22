/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

// Generate elliprd for the JOSS article

use ellip::elliprd;
use ellip_plot_graph::*;
use plotly::{
    color::NamedColor,
    common::{Font, Line, Mode, Title},
    layout::Axis, Layout, Plot, Scatter,
};

macro_rules! get_trace {
    ($x: expr, $y: expr, $name: expr) => {{
        let value = $x
            .iter()
            .map(|&xi| match elliprd(xi, $y, 1.0) {
                Ok(ans) => ans,
                Err(_) => f64::NAN,
            })
            .collect();

        Scatter::new($x.clone(), value)
            .mode(Mode::Lines)
            .name($name)
    }};
    ($x: expr, $y: expr, $name: expr, $line_color: expr) => {
        get_trace!($x, $y, $name).line(Line::new().color($line_color).width(2.0))
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
        get_trace!(&x, 1.0, "y=1.0"),
        get_trace!(&x, 5.0, "y=5.0", NamedColor::Blue),
        get_trace!(&x, 25.0, "y=25.0"),
    ]);
    plot.set_layout(
        Layout::new()
            .x_axis(
                Axis::new()
                    .title(Title::with_text("x").font(Font::new().size(18)))
                    .tick_font(Font::new().size(16))
                    .show_line(true),
            )
            .y_axis(
                Axis::new()
                    .title(Title::with_text("RD(x,y,1)").font(Font::new().size(18)))
                    .tick_font(Font::new().size(16))
                    .show_line(true)
                    .range(vec![0.0, 2.0]),
            )
            .legend(
                plotly::layout::Legend::new()
                    .x(1.0)
                    .x_anchor(plotly::common::Anchor::Right),
            ),
    );

    make_html!(plot, "elliprd_plot.html");
    write_svg!(plot, "elliprd_plot_paper.svg", 1100, 600, 1.0);
    println!("Done");
    Ok(())
}
