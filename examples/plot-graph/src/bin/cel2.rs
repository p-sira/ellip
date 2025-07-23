/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::cel2;
use ellip_plot_graph::*;
use plotly::{
    ImageFormat, Layout, Plot, Scatter,
    color::NamedColor,
    common::{Line, Mode},
    layout::{Annotation, Axis},
};

macro_rules! get_trace {
    ($kc: expr, $a: expr, $b: expr, $name: expr) => {{
        let value = $kc
            .iter()
            .map(|&kci| match cel2(kci, $a, $b) {
                Ok(ans) => ans,
                Err(_) => f64::NAN,
            })
            .collect();

        Scatter::new($kc.clone(), value)
            .mode(Mode::Lines)
            .name($name)
    }};
    ($kc: expr, $a: expr, $b: expr, $name: expr, $line_color: expr) => {
        get_trace!($kc, $a, $b, $name).line(Line::new().color($line_color))
    };
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n_points = 100;
    let range_kc = [-3, 3];
    let mut kc: Vec<f64> = [
        (range_kc[0] * n_points..=range_kc[1] * n_points)
            .map(|x| x as f64 / n_points as f64)
            .collect::<Vec<f64>>(),
        // Make the plot more dense near zero to improve the visual
        (-10..10)
            .map(|x| x as f64 / 1000.0)
            .skip(1)
            .collect::<Vec<f64>>(),
    ]
    .concat();
    kc.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut plot = Plot::new();
    plot.add_traces(vec![
        get_trace!(&kc, 3.0, 3.0, "a=3, b=3", NamedColor::Red),
        get_trace!(&kc, 1.0, 3.0, "a=1, b=3"),
        get_trace!(&kc, 3.0, 1.0, "a=3, b=1"),
        get_trace!(&kc, 1.0, 1.0, "a=1, b=1", NamedColor::Blue),
    ]);
    plot.set_layout(
        Layout::new()
            .title("Bulirsch's Complete Elliptic Integral of the Second Kind (cel2)")
            .x_axis(Axis::new().title("kc").show_line(true))
            .y_axis(
                Axis::new()
                    .title("cel2(kc,a,b)")
                    .show_line(true)
                    // .range(vec![0.0, 4.0]),
            )
            .annotations(vec![Annotation::new()
            .text(format!(
                "Generated using the function <a href=\"https://docs.rs/ellip/latest/ellip/bulirsch/fn.cel2.html\" target=\"_blank\">cel2</a> from <a href=\"https://crates.io/crates/ellip\" target=\"_blank\">ellip</a> v{}",
                env!("CARGO_PKG_VERSION")
            ))
                .x_ref("paper")
                .y_ref("paper")
                .y(-0.15)
                .x(1.08)
                .show_arrow(false)]),
    );

    make_html!(plot, "cel2_plot.html");
    plot.write_image(
        figure_path!("cel2_plot.svg"),
        ImageFormat::SVG,
        900,
        600,
        1.0,
    );
    Ok(())
}
