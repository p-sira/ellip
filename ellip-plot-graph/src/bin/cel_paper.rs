/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::cel;
use ellip_plot_graph::*;
use plotly::{
    color::NamedColor,
    common::{Font, Line, Mode, Title},
    layout::Axis,
    ImageFormat, Layout, Plot, Scatter,
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
    let range_kc = [-2, 2];
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
            .x_axis(
                Axis::new()
                    .title(Title::with_text("kc").font(Font::new().size(18)))
                    .show_line(true),
            )
            .y_axis(
                Axis::new()
                    .title(Title::with_text("cel(kc,p,1,1)").font(Font::new().size(18)))
                    .show_line(true), // .range(vec![0.0, 4.0]),
            )
            .legend(
                plotly::layout::Legend::new()
                    .x(1.0)
                    .x_anchor(plotly::common::Anchor::Right),
            ),
    );

    make_html!(plot, "cel_plot_paper.html");
    write_svg!(plot, "cel_plot_paper.svg", 800, 600, 1.0);
    println!("Done");
    Ok(())
}
