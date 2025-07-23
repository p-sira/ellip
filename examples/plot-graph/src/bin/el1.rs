/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use std::f64::consts::{FRAC_PI_2, FRAC_PI_3, FRAC_PI_4};

use ellip::el1;
use ellip_plot_graph::*;
use plotly::{
    ImageFormat, Layout, Plot, Scatter,
    color::NamedColor,
    common::{Anchor, Line, Mode},
    layout::{Annotation, Axis, Legend},
};

macro_rules! get_trace {
    ($x: expr, $kc: expr, $name: expr) => {{
        let value = $kc
            .iter()
            .map(|&kci| match el1($x, kci) {
                Ok(ans) => ans,
                Err(_) => f64::NAN,
            })
            .collect();

        Scatter::new($kc.clone(), value)
            .mode(Mode::Lines)
            .name($name)
    }};
    ($x: expr, $kc: expr, $name: expr, $line_color: expr) => {
        get_trace!($x, $kc, $name).line(Line::new().color($line_color))
    };
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n_points = 100;
    let range_kc = [-3, 3];

    let mut kc: Vec<f64> = [
        (range_kc[0] * n_points..=range_kc[1] * n_points)
            .map(|x| x as f64 / n_points as f64)
            .collect(),
        // Make the plot more dense near zero to improve the visual
        (-5..5)
            .map(|x| x as f64 / 500.0)
            .skip(1)
            .collect::<Vec<f64>>(),
    ]
    .concat();
    kc.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut plot = Plot::new();
    plot.add_traces(vec![
        get_trace!(FRAC_PI_2.tan(), &kc, "x=tan(π/2)", NamedColor::Red),
        get_trace!(FRAC_PI_3.tan(), &kc, "x=tan(π/3)"),
        get_trace!(FRAC_PI_4.tan(), &kc, "x=tan(π/4)"),
    ]);

    plot.set_layout(
         Layout::new()
             .title("Bulirsch's Incomplete Elliptic Integral of the First Kind (el1)")
             .x_axis(Axis::new().title("kc").show_line(true))
             .y_axis(
                 Axis::new()
                     .title("el1(x,kc)")
                     .show_line(true)
                     .range(vec![0.0, 5.0]),
             ).legend(Legend::new().y_anchor(Anchor::Middle).y(0.5))
             .annotations(vec![Annotation::new()
             .text(format!(
                 "Generated using the function <a href=\"https://docs.rs/ellip/latest/ellip/bulirsch/fn.el1.html\" target=\"_blank\">el1</a> from <a href=\"https://crates.io/crates/ellip\" target=\"_blank\">ellip</a> v{}",
                 env!("CARGO_PKG_VERSION")
             ))
                 .x_ref("paper")
                 .y_ref("paper")
                 .y(-0.15)
                 .x(1.08)
                 .show_arrow(false)]),
     );

    make_html!("el1_plot.html");
    plot.write_image(
        figure_path!("el1_plot.svg"),
        ImageFormat::SVG,
        900,
        600,
        1.0,
    );
    Ok(())
}
