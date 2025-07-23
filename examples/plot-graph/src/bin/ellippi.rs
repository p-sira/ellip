/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::ellippi;
use ellip_plot_graph::*;
use plotly::{
    Layout, Plot, Surface,
    common::{ColorScale, ColorScalePalette},
    layout::{Annotation, Axis, LayoutScene},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n_points_n = 50;
    let n_points_n_cauchy = 30;
    let n_points_m = 50;
    let range_n = [-2.0, 1.0 - 1e-5];
    let range_n_cauchy = [1.0 + 1e-5, 2.0];
    let range_m = [-2.0, 1.0 - 5.0 * f64::EPSILON];

    let n: Vec<f64> = (0..n_points_n)
        .map(|i| range_n[0] + i as f64 * (range_n[1] - range_n[0]) / (n_points_n - 1) as f64)
        .collect();
    let n_cauchy: Vec<f64> = (0..n_points_n_cauchy)
        .map(|i| {
            range_n_cauchy[0]
                + i as f64 * (range_n_cauchy[1] - range_n_cauchy[0])
                    / (n_points_n_cauchy - 1) as f64
        })
        .collect();
    let m: Vec<f64> = (0..n_points_m)
        .map(|i| range_m[0] + i as f64 * (range_m[1] - range_m[0]) / (n_points_m - 1) as f64)
        .collect();

    let ellippi_values = m
        .iter()
        .map(|&mi| {
            n.iter()
                .map(|&nj| match ellippi(nj, mi) {
                    Ok(ans) => ans,
                    Err(_) => f64::NAN,
                })
                .collect()
        })
        .collect();

    let trace = Surface::new(ellippi_values)
        .x(n)
        .y(m.clone())
        .name("Π(n,m)")
        .color_scale(ColorScale::Palette(ColorScalePalette::Viridis))
        .cmin(-6.0)
        .cmax(6.0);

    let ellippi_cauchy_values = m
        .iter()
        .map(|&mi| {
            n_cauchy
                .iter()
                .map(|&nj| match ellippi(nj, mi) {
                    Ok(ans) => ans,
                    Err(_) => f64::NAN,
                })
                .collect()
        })
        .collect();

    let trace_cauchy = Surface::new(ellippi_cauchy_values)
        .x(n_cauchy)
        .y(m)
        .name("Π(n,m)")
        .color_scale(ColorScale::Palette(ColorScalePalette::Viridis))
        .cmin(-6.0)
        .cmax(6.0);

    let mut plot = Plot::new();
    plot.add_trace(trace);
    plot.add_trace(trace_cauchy);
    plot.set_layout(
        Layout::new()
            .title("Complete Elliptic Integral of the Third Kind (Π)")
            .width(800)
            .height(600)
            .scene(
                LayoutScene::new()
                    .x_axis(Axis::new().title("n").show_line(true))
                    .y_axis(
                        Axis::new()
                            .title("m")
                            .show_line(true)
                            .range(vec![range_m[1], range_m[0]]),
                    )
                    .z_axis(
                        Axis::new()
                            .title("Π(n,m)")
                            .show_line(true)
                            .range(vec![-6.0, 6.0]),
                    )
            )
            .annotations(vec![Annotation::new()
            .text(format!(
                "Generated using the function <a href=\"https://docs.rs/ellip/latest/ellip/legendre/fn.ellippi.html\" target=\"_blank\">ellippi</a> from <a href=\"https://crates.io/crates/ellip\" target=\"_blank\">ellip</a> v{}",
                env!("CARGO_PKG_VERSION")
            ))
            .x_ref("paper")
            .y_ref("paper")
            .y(-0.15)
            .x(1.08)
            .show_arrow(false)]),
    );

    make_html!("ellippi_plot.html");
    // Current plotly.rs doesn't support exporting 3D plot as image.
    // The workaround is using the capture function in the html to save a png file.
    // plot.write_image(figure_path!("ellippi_plot.svg"), ImageFormat::SVG, 900, 900, 0.2);
    Ok(())
}
