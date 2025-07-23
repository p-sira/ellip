/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

macro_rules! check {
    ($fn_name:ident, $check_method:ident, $value_name:expr, [$($var:ident),* $(,)?] $(,)?) => {
        $(
            if $var.$check_method() {
                return Err(concat![stringify!($fn_name), ": ", stringify!($var), " cannot be ", $value_name, "."])
            }
        )*
    };
    ($fn_name:ident, ($predicate:expr), $value_name:expr, [$($var:ident),* $(,)?] $(,)?) => {
        $(
            if $predicate {
                return Err(concat![stringify!($fn_name), ": ", stringify!($var), " cannot be ", $value_name, "."])
            }
        )*
    };
    (@nan, $fn_name:ident, [$($var:ident),* $(,)?] $(,)?) => {
        check!($fn_name, is_nan, "nan", [$($var),*])
    };
    (@zero, $fn_name:ident, [$($var:ident),* $(,)?] $(,)?) => {
        check!($fn_name, is_zero, "zero", [$($var),*])
    }
}
pub(crate) use check;

macro_rules! case {
    ($predicate:expr, $res:expr) => {
        if $predicate {
            return Ok($res);
        }
    };
    (@any [$($var:expr),+] == $value:expr, $res:expr) =>{
        $(
            if $var == $value {
                return Ok($res)
            }
        )+
    }
}
pub(crate) use case;

/// ans is considered valid when ans is finite and not nan.
macro_rules! return_if_valid_else {
    ($ans:expr, {$($else:tt)+}) => {
        if $ans.is_finite() && !$ans.is_nan() {
            return Ok($ans);
        } else {
            $($else)+
            unreachable!()
        }
    };
}
pub(crate) use return_if_valid_else;

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
