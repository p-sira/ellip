/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

#[macro_export]
macro_rules! figure_path {
    ($filename:literal) => {
        #[cfg(not(feature = "save-to-bin"))]
        concat!["../../figures/", $filename]
        #[cfg(feature = "save-to-bin")]
        concat!["../../figures/bin/", $filename]
    };
}

#[macro_export]
macro_rules! make_html {
    ($plot:ident, $filename:literal) => {
        println!("Making HTML for {}", $filename);
        #[cfg(feature = "open-html")]
        $plot.show_html(figure_path!($filename));
        #[cfg(not(feature = "open-html"))]
        $plot.write_html(figure_path!($filename));
    };
}

#[macro_export]
macro_rules! write_svg {
    ($plot:ident, $filename:literal, $width:expr, $height:expr, $scale:expr) => {
        println!("Writing image to {}", $filename);
        $plot.write_image(
            figure_path!($filename),
            ImageFormat::SVG,
            $width,
            $height,
            $scale,
        );
    };
}

pub fn ellip_version() -> String {
    use std::process::Command;

    let output = Command::new("cargo")
        .args(["tree", "--invert", "--package", "ellip"])
        .output()
        .unwrap()
        .stdout;

    String::from_utf8_lossy(&output)
        .lines()
        .next()
        .and_then(|line| line.strip_prefix("ellip v"))
        .unwrap()
        .split_whitespace()
        .next()
        .unwrap()
        .to_owned()
}
