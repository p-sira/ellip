#![cfg(not(feature = "test_force_fail"))]
#![allow(clippy::excessive_precision)]

#[macro_use]
mod common;

macro_rules! assert_close {
    ($expected:expr, $actual:expr, $rtol:expr) => {
        ellip::util::assert_close($actual, $expected, $rtol);
    };
}

mod accuracy;
mod identities;
mod optimization_regression;
