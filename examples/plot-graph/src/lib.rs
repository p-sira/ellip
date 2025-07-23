/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

#[macro_export]
macro_rules! figure_path {
    ($filename:literal) => {
        concat!["../../figures/", $filename]
    };
}

#[macro_export]
macro_rules! make_html {
    ($filename:literal) => {
        #[cfg(feature="open-html")]
        plot.show_html(figure_path!($filename));
        #[cfg(not(feature="open-html"))]
        plot.write_html(figure_path!($filename));
    };
}