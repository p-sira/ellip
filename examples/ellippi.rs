/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::ellippi;
use plotters::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let root = SVGBackend::new("figures/ellippi_plot.svg", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    // Define ranges
    let n_points = 50;
    let range_n = [-2, 1];
    let range_n_over_one = [1, 2];
    let range_m = [-2, 1];
    let range_y = [-6.0, 6.0];

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Complete Elliptic Integral of the Third Kind",
            ("serif", 30),
        )
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_3d(
            range_n_over_one[1] as f64..range_n[0] as f64,
            range_y[0]..range_y[1],
            range_m[1] as f64..range_m[0] as f64,
        )?;

    chart.with_projection(|mut p| {
        p.yaw = 0.14;
        p.into_matrix()
    });

    chart.configure_axes().draw()?;

    // Plot the result
    // Separate the drawing process because the function is discontinuous
    chart.draw_series(
        SurfaceSeries::xoz(
            (range_n[0] * n_points..=range_n[1] * n_points).map(|x| {
                match x as f64 / n_points as f64 {
                    1.0 => 1.0 - f64::EPSILON,
                    i => i,
                }
            }),
            (range_m[0] * n_points..=range_m[1] * n_points).map(|z| {
                match z as f64 / n_points as f64 {
                    1.0 => 1.0 - f64::EPSILON,
                    i => i,
                }
            }),
            |n, m| {
                ellippi(n, m)
                    .unwrap()
                    .clamp(range_y[0] - 0.5, range_y[1] + 0.5)
            },
        )
        .style_func(&|&y| ViridisRGB::get_color_normalized(y, range_y[0], range_y[1]).filled()),
    )?;

    chart.draw_series(
        SurfaceSeries::xoz(
            (range_n_over_one[0] * n_points..=range_n_over_one[1] * n_points).map(|x| {
                match x as f64 / n_points as f64 {
                    1.0 => 1.0 + f64::EPSILON,
                    n => n,
                }
            }),
            (range_m[0] * n_points..=range_m[1] * n_points).map(|z| {
                match z as f64 / n_points as f64 {
                    1.0 => 1.0 - f64::EPSILON,
                    i => i,
                }
            }),
            |n, m| {
                ellippi(n, m)
                    .unwrap()
                    .clamp(range_y[0] - 0.5, range_y[1] + 0.5)
            },
        )
        .style_func(&|&y| ViridisRGB::get_color_normalized(y, range_y[0], range_y[1]).filled()),
    )?;

    root.present()?;
    Ok(())
}
