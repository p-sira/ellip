/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use std::f64::consts::{FRAC_PI_2, FRAC_PI_3, FRAC_PI_4};

use ellip::ellipeinc;
use plotters::{prelude::*, style::full_palette::ORANGE};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let root = SVGBackend::new("figures/ellipeinc_plot.svg", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Incomplete Elliptic Integral of the Second Kind",
            ("serif", 30),
        )
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(-5.0..2.0, 0.0..4.0)?;

    chart
        .configure_mesh()
        .x_desc("m (k²)")
        .y_desc("ellipeinc(φ,m)")
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
            match ellipeinc(phi, m) {
                Ok(res) => points.push((m, res)),
                Err(_) => (),
            }
        }
        points
    }

    let ellipeincpi_2_points: Vec<(f64, f64)> = get_points(FRAC_PI_2);
    let ellipeincpi_3_points: Vec<(f64, f64)> = get_points(FRAC_PI_3);
    let ellipeincpi_4_points: Vec<(f64, f64)> = get_points(FRAC_PI_4);

    // Plot the result
    chart
        .draw_series(LineSeries::new(ellipeincpi_2_points, RED.stroke_width(2)))?
        .label("φ=π/2")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED.stroke_width(2)));

    chart
        .draw_series(LineSeries::new(
            ellipeincpi_3_points,
            ORANGE.stroke_width(2),
        ))?
        .label("φ=π/3")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], ORANGE.stroke_width(2)));

    chart
        .draw_series(LineSeries::new(ellipeincpi_4_points, GREEN.stroke_width(2)))?
        .label("φ=π/4")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], GREEN.stroke_width(2)));

    chart
        .configure_series_labels()
        .border_style(&BLACK)
        .background_style(&WHITE.mix(0.8))
        .label_font(("serif", 25).into_font())
        .position(SeriesLabelPosition::UpperRight)
        .draw()?;

    root.present()?;
    Ok(())
}
