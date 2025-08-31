/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::ellipeinc;
use ellip_plot_graph::*;
use plotly::{
    Layout, Plot, Surface,
    common::{ColorScale, ColorScalePalette},
    layout::{Annotation, AspectRatio, Axis, LayoutScene},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n_points_m = 50;
    let n_points_s2p = 50;
    let range_m = [-1.0, 1.0];
    let range_s2p = [0.0, 0.5];

    let mut m: Vec<f64> = (0..n_points_m)
        .map(|i| range_m[0] + i as f64 * (range_m[1] - range_m[0]) / (n_points_m - 1) as f64)
        .collect();

    let mut s2p: Vec<f64> = (0..n_points_s2p)
        .map(|i| {
            range_s2p[0] + i as f64 * (range_s2p[1] - range_s2p[0]) / (n_points_s2p - 1) as f64
        })
        .collect();

    let n_points_dense = 500;
    m.extend(
        (0..n_points_dense)
            .map(|i| 1.0 + i as f64 * (2.0 - 1.0) / (n_points_dense - 1) as f64)
            .collect::<Vec<f64>>(),
    );
    s2p.extend(
        (0..n_points_dense)
            .map(|i| 0.5 + i as f64 * (1.0 - 0.5) / (n_points_dense - 1) as f64)
            .collect::<Vec<f64>>(),
    );

    let ellipeinc_values: Vec<Vec<f64>> = s2p
        .iter()
        .map(|&s2pi| {
            m.iter()
                .map(|&mj| {
                    let phi = s2pi.sqrt().asin();
                    match ellipeinc(phi, mj) {
                        Ok(ans) => ans,
                        Err(_) => f64::NAN,
                    }
                })
                .collect()
        })
        .collect();

    let trace = Surface::new(ellipeinc_values)
        .x(m)
        .y(s2p)
        .name("E(φ,m)")
        .color_scale(ColorScale::Palette(ColorScalePalette::Viridis))
        .cmin(0.0)
        .cmax(2.0);

    let mut plot = Plot::new();
    plot.add_trace(trace);
    plot.set_layout(
        Layout::new()
            .title("Incomplete Elliptic Integral of the Second Kind (E)")
            .width(800)
            .height(600)
            .scene(
                LayoutScene::new()
                    .x_axis(Axis::new().title("m").show_line(true))
                    .y_axis(
                        Axis::new()
                            .title("sin²φ")
                            .show_line(true)
                    )
                    .z_axis(
                        Axis::new()
                            .title("E(φ,m)")
                            .show_line(true)
                            .range(vec![0.0, 2.5]),
                    ).aspect_ratio(AspectRatio::new().x(1.0).y(1.0))
            )
            .annotations(vec![Annotation::new()
            .text(format!(
                "Generated using <a href=\"https://docs.rs/ellip/latest/ellip/legendre/fn.ellipeinc.html\" target=\"_blank\">ellipeinc</a> from <a href=\"https://crates.io/crates/ellip\" target=\"_blank\">ellip</a> v{}",
                ellip_version()
            ))
            .x_ref("paper")
            .y_ref("paper")
            .y(-0.15)
            .x(1.08)
            .show_arrow(false)]),
    );

    make_html!(plot, "ellipeinc_plot_3d.html");
    // Current plotly.rs doesn't support exporting 3D plot as image.
    // The workaround is using the capture function in the html to save a png file.
    // plot.write_image(figure_path!("ellipeinc_plot_3d.svg"), ImageFormat::SVG, 900, 900, 0.2);
    Ok(())
}
