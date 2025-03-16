/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

use ellip::elliprj;
use plotters::{prelude::*, style::full_palette::ORANGE};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let root = SVGBackend::new("figures/elliprj_plot.svg", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    // Compute the integral
    fn get_points(y: f64) -> Vec<(f64, f64)> {
        let n_points = 500;
        let range = [0, 1];
        let mut points: Vec<(f64, f64)> = Vec::with_capacity(n_points);
        for i in range[0] * n_points as i32..range[1] * n_points as i32 {
            let x = i as f64 / n_points as f64;
            match elliprj(x, y, 1.0, 2.0) {
                Ok(res) => points.push((x, res)),
                Err(_) => (),
            }
        }
        points
    }

    let rj_0_1_points = get_points(0.0);
    let rj_01_1_points = get_points(0.1);
    let rj_05_1_points = get_points(0.5);
    let rj_1_1_points = get_points(1.0);

    fn get_neg_points(y: f64) -> Vec<(f64, f64)> {
        let n_points = 500;
        let range = [0, 1];
        let mut points: Vec<(f64, f64)> = Vec::with_capacity(n_points);
        for i in range[0] * n_points as i32..range[1] * n_points as i32 {
            let x = i as f64 / n_points as f64;
            match elliprj(x, y, 1.0, -2.0) {
                Ok(res) => {
                    if res >= -5.0 {
                        points.push((x, res))
                    }
                }
                Err(_) => (),
            }
        }
        points
    }

    let rj_0_1_neg_points = get_neg_points(0.0);
    let rj_01_1_neg_points = get_neg_points(0.1);
    let rj_05_1_neg_points = get_neg_points(0.5);
    let rj_1_1_neg_points = get_neg_points(1.0);

    // Plot the result
    {
        let mut chart = ChartBuilder::on(&root)
            .caption(
                "Symmetric Elliptic Integral of the Third Kind (elliprj)",
                ("serif", 30),
            )
            .margin(20)
            .x_label_area_size(50)
            .y_label_area_size(60)
            .build_cartesian_2d(0.0..1.0, -5.0..5.0)?;

        chart
            .configure_mesh()
            .x_desc("x")
            .y_desc("Value")
            .axis_desc_style(("serif", 25).into_font())
            .label_style(("serif", 20).into_font())
            .draw()?;

        chart
            .draw_series(LineSeries::new(rj_0_1_points, RED.stroke_width(2)))?
            .label("elliprj(x,0,1,2)")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

        chart
            .draw_series(LineSeries::new(rj_01_1_points, ORANGE.stroke_width(2)))?
            .label("elliprj(x,0.1,1,2)")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &ORANGE));

        chart
            .draw_series(LineSeries::new(rj_05_1_points, GREEN.stroke_width(2)))?
            .label("elliprj(x,0.5,1,2)")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

        chart
            .draw_series(LineSeries::new(rj_1_1_points, BLUE.stroke_width(2)))?
            .label("elliprj(x,1,1,2)")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

        chart
            .configure_series_labels()
            .border_style(&BLACK)
            .background_style(&WHITE.mix(0.8))
            .label_font(("serif", 20).into_font())
            .position(SeriesLabelPosition::UpperRight)
            .draw()?;

        // Draw y=0 line for separation
        chart.draw_series(LineSeries::new(
            (0..=10)
                .map(|x| (x as f64 / 10.0, 0.0))
                .collect::<Vec<(f64, f64)>>(),
            BLACK.stroke_width(1),
        ))?;
    }

    // Negative p
    {
        let mut chart = ChartBuilder::on(&root)
            .caption(
                "Symmetric Elliptic Integral of the Third Kind (elliprj)",
                ("serif", 30).with_color(TRANSPARENT),
            )
            .margin(20)
            .x_label_area_size(50)
            .y_label_area_size(60)
            .build_cartesian_2d(0.0..1.0, -5.0..5.0)?;

        chart
            .configure_mesh()
            .x_desc("x")
            .y_desc("Value")
            .axis_desc_style(("serif", 25).into_font().color(&TRANSPARENT))
            .label_style(("serif", 20).into_font().color(&TRANSPARENT))
            .disable_x_axis()
            .disable_x_mesh()
            .disable_y_axis()
            .disable_y_mesh()
            .draw()?;

        chart
            .draw_series(LineSeries::new(rj_0_1_neg_points, RED.stroke_width(2)))?
            .label("elliprj(x,0,1,-2)")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

        chart
            .draw_series(LineSeries::new(rj_01_1_neg_points, ORANGE.stroke_width(2)))?
            .label("elliprj(x,0.1,1,-2)")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &ORANGE));

        chart
            .draw_series(LineSeries::new(rj_05_1_neg_points, GREEN.stroke_width(2)))?
            .label("elliprj(x,0.5,1,-2)")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

        chart
            .draw_series(LineSeries::new(rj_1_1_neg_points, BLUE.stroke_width(2)))?
            .label("elliprj(x,1,1,-2)")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

        chart
            .configure_series_labels()
            .border_style(&BLACK)
            .background_style(&WHITE.mix(0.8))
            .label_font(("serif", 20).into_font())
            .position(SeriesLabelPosition::LowerRight)
            .draw()?;
    }

    root.present()?;
    Ok(())
}
