# Testing
This report presents the accuracy of the ellip crate using [symmetric relative error](https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/relative_error.html)
metric. The errors are expressed in units of machine epsilon (ε).
The reference values are computed using [Wolfram Engine](https://www.wolfram.com/engine/),
You can find the scripts in the directory [tests/wolfram/](https://github.com/p-sira/ellip/blob/main/tests/wolfram/).
This report is generated on x86_64-unknown-linux-gnu rustc 1.83.0 using ellip v0.3.1 at `f64` precision (ε≈2.22e-16).

## Legendre's Complete Elliptic Integrals
| Function              | Mean (ε) | Median (ε) | Variance (ε²) | Max (ε) |
|-----------------------|----------|------------|---------------|---------|
| ellipk                | 0.93     | 0.62       | 2.09          | 20.09   |
| ellipk (Negative m)   | 1.35     | 0.89       | 34.58         | 122.33  |
| ellipk (Near limits)  | 119.25   | 98.03      | 4019.08       | 270.99  |
| ellipe                | 0.92     | 0.77       | 0.14          | 3.79    |
| ellipe (Negative m)   | 1.05     | 0.89       | 0.22          | 3.12    |
| ellipe (Near limits)  | 1.55     | 1.00       | 0.86          | 6.00    |
| ellippi               | 0.96     | 0.81       | 0.30          | 10.88   |
| ellippi (Negative m)  | 0.91     | 0.81       | 0.14          | 3.17    |
| ellippi (p.v.)        | 1.27     | 0.97       | 0.76          | 13.25   |
| ellippi (Near limits) | 39.90    | 0.97       | 6424.03       | 386.32  |

## Legendre's Incomplete Elliptic Integrals
| Function                | Mean (ε)   | Median (ε) | Variance (ε²)      | Max (ε)      |
|-------------------------|------------|------------|--------------------|--------------|
| ellipf                  | 15.28      | 10.88      | 240.43             | 295.66       |
| ellipf (Negative m)     | 11.73      | 9.76       | 79.05              | 45.10        |
| ellipf (Near limits)    | 1895945.32 | 50.07      | 93450031767571.67  | 51187940.87  |
| ellipeinc               | 10.72      | 8.68       | 72.94              | 53.48        |
| ellipeinc (Negative m)  | 12.74      | 10.23      | 98.71              | 50.54        |
| ellipeinc (Near limits) | 5322751.06 | 49.52      | 736605879806874.38 | 143712609.34 |

## Bulirsch's Complete Elliptic Integrals

## Bulirsch's Incomplete Elliptic Integrals

## Carlson's Symmetric Elliptic Integrals

