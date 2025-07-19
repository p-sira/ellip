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
