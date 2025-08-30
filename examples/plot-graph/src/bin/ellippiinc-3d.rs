/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::ellippiinc;
use ellip_plot_graph::*;
use plotly::{
    common::{ColorScale, ColorScalePalette},
    layout::{Annotation, AspectRatio, Axis, LayoutScene},
    Layout, Plot, Surface,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // s2p < 0.5
    let trace_lt_05;
    {
        let n_points_m = 100;
        let n_points_s2p = 50;
        let range_m = [-1.0, 1.5];
        let range_s2p = [0.0, 0.3];

        let mut m: Vec<f64> = (0..n_points_m)
            .map(|i| range_m[0] + i as f64 * (range_m[1] - range_m[0]) / (n_points_m - 1) as f64)
            .collect();

        let mut s2p: Vec<f64> = (0..n_points_s2p)
            .map(|i| {
                range_s2p[0] + i as f64 * (range_s2p[1] - range_s2p[0]) / (n_points_s2p - 1) as f64
            })
            .collect();

        let n_points_dense = 600;
        m.extend(
            (0..n_points_dense)
                .map(|i| 1.5 + i as f64 * (3.0 - 1.5) / (n_points_dense - 1) as f64)
                .collect::<Vec<f64>>(),
        );
        s2p.extend(
            (0..n_points_dense)
                .map(|i| 0.3 + i as f64 * (0.5 - 0.3) / (n_points_dense - 1) as f64)
                .collect::<Vec<f64>>(),
        );

        let ellippiinc_values: Vec<Vec<f64>> = s2p
            .iter()
            .map(|&s2pi| {
                m.iter()
                    .map(|&mj| {
                        let phi = s2pi.sqrt().asin();
                        match ellippiinc(phi, 2.0, mj) {
                            Ok(ans) => ans,
                            Err(_) => f64::NAN,
                        }
                    })
                    .collect()
            })
            .collect();

        trace_lt_05 = Surface::new(ellippiinc_values)
            .x(m)
            .y(s2p)
            .name("Π(φ,2,m)")
            .color_scale(ColorScale::Palette(ColorScalePalette::Viridis))
            .cmin(0.0)
            .cmax(5.0);
    }

    // s2p > 0.5
    let trace_gt_05;
    {
        let n_points_m = 100;
        let n_points_s2p = 800;
        let range_m = [-1.0, 3.0];
        let range_s2p = [0.5, 1.0];

        let mut m: Vec<f64> = (0..n_points_m)
            .map(|i| range_m[0] + i as f64 * (range_m[1] - range_m[0]) / (n_points_m - 1) as f64)
            .collect();

        let s2p: Vec<f64> = (0..n_points_s2p)
            .map(|i| {
                range_s2p[0] + i as f64 * (range_s2p[1] - range_s2p[0]) / (n_points_s2p - 1) as f64
            })
            .collect();

        let n_points_dense = 700;
        m.extend(
            (0..n_points_dense)
                .map(|i| 0.75 + i as f64 * (2.0 - 0.75) / (n_points_dense - 1) as f64)
                .collect::<Vec<f64>>(),
        );

        let ellippiinc_values: Vec<Vec<f64>> = s2p
            .iter()
            .map(|&s2pi| {
                m.iter()
                    .map(|&mj| {
                        let phi = s2pi.sqrt().asin();
                        match ellippiinc(phi, 2.0, mj) {
                            Ok(ans) => ans,
                            Err(_) => f64::NAN,
                        }
                    })
                    .collect()
            })
            .collect();

        trace_gt_05 = Surface::new(ellippiinc_values)
            .x(m)
            .y(s2p)
            .name("Π(φ,2,m)")
            .color_scale(ColorScale::Palette(ColorScalePalette::Viridis))
            .cmin(0.0)
            .cmax(5.0);
    }

    let mut plot = Plot::new();
    plot.add_trace(trace_lt_05);
    plot.add_trace(trace_gt_05);
    plot.set_layout(
         Layout::new()
             .title("Incomplete Elliptic Integral of the Third Kind (Π)")
             .width(800)
             .height(600)
             .scene(
                 LayoutScene::new()
                     .x_axis(Axis::new().title("m").show_line(true))
                     .y_axis(
                         Axis::new()
                             .title("sin²φ")
                             .show_line(true)
                             .range(vec![0.0, 1.0])
                     )
                     .z_axis(
                         Axis::new()
                             .title("Π(φ,2,m)")
                             .show_line(true)
                             .range(vec![-3.0, 5.0]),
                     ).aspect_ratio(AspectRatio::new().x(1.0).y(1.0))
             )
             .annotations(vec![Annotation::new()
             .text(format!(
                 "Generated using the function <a href=\"https://docs.rs/ellip/latest/ellip/legendre/fn.ellippiinc.html\" target=\"_blank\">ellippiinc</a> from <a href=\"https://crates.io/crates/ellip\" target=\"_blank\">ellip</a> v{}",
                 ellip_version()
             ))
             .x_ref("paper")
             .y_ref("paper")
             .y(-0.15)
             .x(1.08)
             .show_arrow(false)]),
     );

    make_html!(plot, "ellippiinc_plot_3d.html");
    // Current plotly.rs doesn't support exporting 3D plot as image.
    // The workaround is using the capture function in the html to save a png file.
    // plot.write_image(figure_path!("ellippiinc_plot_3d.svg"), ImageFormat::SVG, 900, 900, 0.2);
    Ok(())
}
