/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::elliprf;
use plotters::{prelude::*, style::full_palette::ORANGE};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let root = SVGBackend::new("figures/elliprf_plot.svg", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Symmetric Elliptic Integral of the First Kind (elliprf)",
            ("serif", 30),
        )
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(0.0..1.0, 1.0..5.0)?;

    chart
        .configure_mesh()
        .x_desc("x")
        .y_desc("Value")
        .axis_desc_style(("serif", 25).into_font())
        .label_style(("serif", 20).into_font())
        .draw()?;

    // Compute the integral
    fn get_points(y: f64) -> Vec<(f64, f64)> {
        let n_points = 500;
        let range = [0, 1];
        let mut points: Vec<(f64, f64)> = Vec::with_capacity(n_points);
        for i in range[0] * n_points as i32..range[1] * n_points as i32 {
            let x = i as f64 / n_points as f64;
            match elliprf(x, y, 1.0) {
                Ok(res) => points.push((x, res)),
                Err(_) => (),
            }
        }
        points
    }

    let rf_0_1_points = get_points(0.0);
    let rf_01_1_points = get_points(0.1);
    let rf_05_1_points = get_points(0.5);
    let rf_1_1_points = get_points(1.0);

    // Plot the result
    chart
        .draw_series(LineSeries::new(rf_0_1_points, RED.stroke_width(2)))?
        .label("elliprf(x,0,1)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .draw_series(LineSeries::new(rf_01_1_points, ORANGE.stroke_width(2)))?
        .label("elliprf(x,0.1,1)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &ORANGE));

    chart
        .draw_series(LineSeries::new(rf_05_1_points, GREEN.stroke_width(2)))?
        .label("elliprf(x,0.5,1)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

    chart
        .draw_series(LineSeries::new(rf_1_1_points, BLUE.stroke_width(2)))?
        .label("elliprf(x,1,1)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    chart
        .configure_series_labels()
        .border_style(&BLACK)
        .background_style(&WHITE.mix(0.8))
        .label_font(("serif", 20).into_font())
        .position(SeriesLabelPosition::UpperRight)
        .draw()?;

    root.present()?;
    Ok(())
}
