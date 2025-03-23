/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::elliprc;
use plotters::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let root = SVGBackend::new("figures/elliprc_plot.svg", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Degenerate Symmetric Elliptic Integral of RF (RC)",
            ("serif", 30),
        )
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(-0.0..5.0, 0.0..2.0)?;

    chart
        .configure_mesh()
        .x_desc("x")
        .y_desc("elliprc(x,y)")
        .axis_desc_style(("serif", 25).into_font())
        .label_style(("serif", 20).into_font())
        .draw()?;

    // Compute the integral
    let n_points = 100;
    let range = [0, 5];
    let elliprc1_points: Vec<(f64, f64)> = (range[0] * n_points..range[1] * n_points)
        .map(|i| {
            let x = i as f64 / n_points as f64;
            (x, elliprc(x, 1.0).unwrap())
        })
        .collect();

    let elliprcneg1_points: Vec<(f64, f64)> = (range[0] * n_points..range[1] * n_points)
        .map(|i| {
            let x = i as f64 / n_points as f64;
            (x, elliprc(x, -1.0).unwrap())
        })
        .collect();

    // Plot the result
    chart
        .draw_series(LineSeries::new(elliprc1_points, BLUE.stroke_width(2)))?
        .label("y=1")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE.stroke_width(2)));
    chart
        .draw_series(LineSeries::new(elliprcneg1_points, RED.stroke_width(2)))?
        .label("y=-1")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED.stroke_width(2)));

    chart
        .configure_series_labels()
        .border_style(&BLACK)
        .background_style(&WHITE.mix(0.8))
        .label_font(("serif", 25).into_font())
        .position(SeriesLabelPosition::LowerRight)
        .draw()?;

    root.present()?;
    Ok(())
}
