/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::ellipdinc;
use ellip_plot_graph::*;
use plotly::{
    Layout, Plot, Surface,
    common::{ColorScale, ColorScalePalette},
    layout::{Annotation, AspectRatio, Axis, Camera, CameraCenter, Eye, LayoutScene},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n_points_m = 50;
    let n_points_s2p = 50;
    let range_m = [-1.0, 0.8];
    let range_s2p = [0.0, 0.45];

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
            .map(|i| 0.8 + i as f64 * (2.0 - 0.8) / (n_points_dense - 1) as f64)
            .collect::<Vec<f64>>(),
    );
    s2p.extend(
        (0..n_points_dense)
            .map(|i| 0.45 + i as f64 * (1.0 - 0.45) / (n_points_dense - 1) as f64)
            .collect::<Vec<f64>>(),
    );

    let ellipdinc_values: Vec<Vec<f64>> = s2p
        .iter()
        .map(|&s2pi| {
            m.iter()
                .map(|&mj| {
                    let phi = s2pi.sqrt().asin();
                    match ellipdinc(phi, mj) {
                        Ok(ans) => ans,
                        Err(_) => f64::NAN,
                    }
                })
                .collect()
        })
        .collect();

    let trace = Surface::new(ellipdinc_values)
        .x(m)
        .y(s2p)
        .name("D(φ,m)")
        .color_scale(ColorScale::Palette(ColorScalePalette::Viridis))
        .cmin(0.0)
        .cmax(3.5);

    let mut plot = Plot::new();
    plot.add_trace(trace);
    plot.set_layout(
        Layout::new()
            .title("Incomplete Elliptic Integral of the Legendre's Type (D)")
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
                            .title("D(φ,m)")
                            .show_line(true)
                            .range(vec![0.0, 3.5]),
                    )
                    .aspect_ratio(AspectRatio::new().x(1.0).y(1.0))
                    .camera(
                        Camera::new().center(
                            CameraCenter::from((-0.004331696578336765, -0.05157520610852415, -0.21574139614464038))
                        )
                        .eye(Eye::from((1.430571919533991, -1.6290477053961274, 0.15860086230220793)))
                    )
            )
            .legend(
                plotly::layout::Legend::new()
                    .x(1.0)
                    .x_anchor(plotly::common::Anchor::Right)
            )
            .annotations(vec![Annotation::new()
            .text(format!(
                "Generated using <a href=\"https://docs.rs/ellip/latest/ellip/legendre/fn.ellipdinc.html\" target=\"_blank\">ellipdinc</a> from <a href=\"https://crates.io/crates/ellip\" target=\"_blank\">ellip</a> v{}",
                ellip_version()
            ))
            .x_ref("paper")
            .y_ref("paper")
            .y(-0.15)
            .x(1.08)
            .show_arrow(false)]),
    );

    make_html!(plot, "ellipdinc_3d.html");
    write_svg!(plot, "ellipdinc_3d.svg", 800, 600, 1.0);
    Ok(())
}
