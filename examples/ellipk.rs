/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::ellipk;
use plotters::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let root = SVGBackend::new("figures/ellipk_plot.svg", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Complete Elliptic Integral of the First Kind (ellipk)",
            ("serif", 30),
        )
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(-2.0..1.0, 0.0..4.0)?;

    chart
        .configure_mesh()
        .x_desc("m (k²)")
        .y_desc("ellipk(m)")
        .axis_desc_style(("serif", 25).into_font())
        .label_style(("serif", 20).into_font())
        .draw()?;

    // Compute the integral
    let n_points = 100;
    let range = [-2, 1];
    let ellipk_points: Vec<(f64, f64)> = (range[0] * n_points..range[1] * n_points)
        .map(|x| {
            let m = x as f64 / n_points as f64;
            (m, ellipk(m).unwrap())
        })
        .collect();

    // Plot the result
    chart.draw_series(LineSeries::new(ellipk_points, RED.stroke_width(2)))?;
    
    root.present()?;
    Ok(())
}
