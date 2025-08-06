/*
 * Ellip is licensed under The 3-Clause BSD, see LICENSE.
 * Copyright 2025 Sira Pornsiriprasert <code@psira.me>
 */

/// Macro to conditionally return error.
macro_rules! check {
    (@return Err, $fn_name:ident, $var:ident, $value_name:expr) => {
        return Err(concat![stringify!($fn_name), ": ", stringify!($var), " cannot be ", $value_name, "."])
    };
    ($fn_name:ident, $check_method:ident, $value_name:expr, [$($var:ident),* $(,)?] $(,)?) => {
        $(
            if $var.$check_method() {
                check!(@return Err, $fn_name, $var, $value_name);
            }
        )*
    };
    ($fn_name:ident, ($predicate:expr), $value_name:expr, [$($var:ident),* $(,)?] $(,)?) => {
        $(
            if $predicate {
                check!(@return Err, $fn_name, $var, $value_name);
            }
        )*
    };
    (@nan, $fn_name:ident, [$($var:ident),* $(,)?] $(,)?) => {{
        $(
            if $var.is_nan() {
                return Err(concat![stringify!($fn_name), ": Arguments cannot be NAN."])
            }
        )*

    }};
    (@zero, $fn_name:ident, [$($var:ident),* $(,)?] $(,)?) => {
        check!($fn_name, is_zero, "zero", [$($var),*])
    };
    (@inf, $fn_name:ident, [$($var:ident),* $(,)?] $(,)?) => {
        check!($fn_name, is_infinite, "infinite", [$($var),*])
    };
    (@neg, $fn_name:ident, $msg:expr, [$first:ident, $($var:ident),* $(,)?] $(,)?) => {
        if $first$(.min($var))* < T::zero() {
            return Err(concat![stringify!($fn_name), ": ", $msg]);
        }
    };
    (@neg, $fn_name:ident, [$first:ident, $($var:ident),* $(,)?] $(,)?) => {
        check!(@neg, $fn_name, "Arguments must be non-negative.", [$first, $($var),*])
    };
    (@multi_zero, $fn_name:ident, [$($var:ident),* $(,)?] $(,)?) => {{
        let mut count: u8 = 0;
        $(
            count += $var.is_zero() as u8;
        )*
        if count > 1 {
            return Err(concat![stringify!($fn_name), ": At most one argument can be zero."]);
        }
    }};
    (@multi, $fn_name:ident, $value_name:expr, $check_method:ident, [$first:ident, $($var:ident),* $(,)?] $(,)?) => {{
        // Check if more than one var leads to true and return error.
        let var_array = [$first, $($var),*];
        let mut count = 0;
        for var in var_array {
            if var.$check_method() {
                count += 1;
            }
            if count >= 2 {
                return Err(
                    concat![
                        stringify!($fn_name),
                        ": More than one argument in ",
                        stringify!($first), $(", ", stringify!($var), )*
                        " cannot be ",
                        $value_name,
                        "."
                    ]
                )
            }
        }
    }};
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

macro_rules! declare {
    (mut [$($var:ident),+ $(,)?] = $value:expr) => {
        $(let mut $var = $value;)+
    };
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
