/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

macro_rules! check_nan {
    ($fn_name:ident, [$first:ident]) => {
        if $first.is_nan() {
            return Err(concat![stringify!($fn_name), ": Arguments cannot be nan."])
        }
    };
    ($fn_name:ident, [$first:ident, $($other:ident),* $(,)?] $(,)?) => {
        if $first.is_nan() $(|| $other.is_nan())* {
            return Err(concat![stringify!($fn_name), ": Arguments cannot be nan."])
        }
    };
}
pub(crate) use check_nan;

macro_rules! declare {
    (mut [$($var:ident $(= $value:expr)?),+ $(,)?]) => {
        $(let mut $var $(= $value)?;)+
    };
    ([$($var:ident $(= $value:expr)?),+ $(,)?]) => {
        $(let $var $(= $value)?;)+
    };
}
pub(crate) use declare;

macro_rules! let_mut {
    ($($var:ident),+ $(,)?) => {
        $(let mut $var = $var;)+
    };
}
pub(crate) use let_mut;
