/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::cel;
use ellip_plot_graph::*;
use plotly::{
    ImageFormat, Layout, Plot, Scatter,
    color::NamedColor,
    common::{Line, Mode},
    layout::{Annotation, Axis},
};

macro_rules! get_trace {
    ($kc: expr, $p: expr, $name: expr) => {{
        let value = $kc
            .iter()
            .map(|&kci| match cel(kci, $p, 1.0, 1.0) {
                Ok(ans) => ans,
                Err(_) => f64::NAN,
            })
            .collect();

        Scatter::new($kc.clone(), value)
            .mode(Mode::Lines)
            .name($name)
    }};
    ($kc: expr, $p: expr, $name: expr, $line_color: expr) => {
        get_trace!($kc, $p, $name).line(Line::new().color($line_color))
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
        get_trace!(&kc, 1.0, "p=1", NamedColor::Red),
        get_trace!(&kc, 0.5, "p=0.5"),
        get_trace!(&kc, 0.25, "p=0.25"),
        get_trace!(&kc, -1.0, "p=-1", NamedColor::Purple),
        get_trace!(&kc, -0.5, "p=-0.5", NamedColor::DarkViolet),
        get_trace!(&kc, -0.25, "p=-0.25", NamedColor::Blue),
    ]);
    plot.set_layout(
        Layout::new()
            .title("General Complete Elliptic Integral (cel)")
            .x_axis(Axis::new().title("kc").show_line(true))
            .y_axis(
                Axis::new()
                    .title("cel(kc,p,1,1)")
                    .show_line(true)
                    // .range(vec![0.0, 4.0]),
            )
            .annotations(vec![Annotation::new()
            .text(format!(
                "Generated using the function <a href=\"https://docs.rs/ellip/latest/ellip/bulirsch/fn.cel.html\" target=\"_blank\">cel</a> from <a href=\"https://crates.io/crates/ellip\" target=\"_blank\">ellip</a> v{}",
                ellip_version()
            ))
                .x_ref("paper")
                .y_ref("paper")
                .y(-0.15)
                .x(1.08)
                .show_arrow(false)]),
    );

    make_html!(plot, "cel_plot.html");
    write_svg!(plot, "cel_plot.svg", 900, 600, 1.0);
    println!("Done");
    Ok(())
}
