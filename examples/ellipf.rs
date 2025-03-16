/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use std::f64::consts::{FRAC_PI_2, FRAC_PI_3, FRAC_PI_4};

use ellip::ellipf;
use plotters::{prelude::*, style::full_palette::ORANGE};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let root = SVGBackend::new("figures/ellipf_plot.svg", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Incomplete Elliptic Integral of the First Kind (ellipf)",
            ("serif", 30),
        )
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(-5.0..2.0, 0.0..4.0)?;

    chart
        .configure_mesh()
        .x_desc("m (k²)")
        .y_desc("Value")
        .axis_desc_style(("serif", 25).into_font())
        .label_style(("serif", 20).into_font())
        .draw()?;

    // Compute the integral
    fn get_points(phi: f64) -> Vec<(f64, f64)> {
        let n_points = 100;
        let range = [-5, 2];
        let mut points: Vec<(f64, f64)> = Vec::with_capacity(n_points);
        for i in range[0] * n_points as i32..range[1] * n_points as i32 {
            let m = i as f64 / n_points as f64;
            match ellipf(phi, m) {
                Ok(res) => points.push((m, res)),
                Err(_) => (),
            }
        }
        points
    }

    let ellipfpi_2_points: Vec<(f64, f64)> = get_points(FRAC_PI_2);
    let ellipfpi_3_points: Vec<(f64, f64)> = get_points(FRAC_PI_3);
    let ellipfpi_4_points: Vec<(f64, f64)> = get_points(FRAC_PI_4);

    // Plot the result
    chart
        .draw_series(LineSeries::new(ellipfpi_2_points, RED.stroke_width(2)))?
        .label("ellipf(π/2,m)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .draw_series(LineSeries::new(ellipfpi_3_points, ORANGE.stroke_width(2)))?
        .label("ellipf(π/3,m)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &ORANGE));

    chart
        .draw_series(LineSeries::new(ellipfpi_4_points, GREEN.stroke_width(2)))?
        .label("ellipf(π/4,m)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

    chart
        .configure_series_labels()
        .border_style(&BLACK)
        .background_style(&WHITE.mix(0.8))
        .label_font(("serif", 20).into_font()).position(SeriesLabelPosition::LowerRight)
        .draw()?;

    root.present()?;
    Ok(())
}
