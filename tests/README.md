# Testing
This report presents the accuracy of the ellip crate using [**symmetric relative error**](https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/relative_error.html)
metric. Errors are expressed in units of machine epsilon (ε).
The test data spans the domain of each function up to **μ** to avoid approaching the function's limit.
The reference values are computed using [**Wolfram Engine**](https://www.wolfram.com/engine/).
You can find the scripts in the directory [tests/wolfram/](https://github.com/p-sira/ellip/blob/main/tests/wolfram/).
This report is generated on x86_64-unknown-linux-gnu rustc 1.88.0 using ellip v0.3.6 at `f64` precision (ε≈2.22e-16).

## Legendre's Complete Elliptic Integrals
| Function        | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|-----------------|----------|------------|---------|---------|---------------|--------|
| ellipk          | 0.51     | 0.00       | 1.88    | 131.90  | 18.15         | 1      |
| ellipk (Neg m)  | 0.70     | 0.52       | 2.20    | 107.50  | 15.78         | 1      |
| ellipe          | 0.55     | 0.67       | 2.29    | 3.00    | 0.30          | 1      |
| ellipe (Neg m)  | 0.66     | 0.75       | 2.24    | 3.03    | 0.35          | 1      |
| ellippi         | 24.39    | 0.00       | 881.35  | 1216.37 | 16052.14      | 50     |
| ellippi (Neg m) | 65.04    | 0.55       | 2470.49 | 3197.93 | 1.30e5        | 50     |
| ellippi (p.v.)  | 3.36     | 0.71       | 45.99   | 520.23  | 582.48        | 50     |

## Legendre's Incomplete Elliptic Integrals
| Function           | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|--------------------|----------|------------|---------|---------|---------------|--------|
| ellipf             | 0.81     | 0.72       | 3.27    | 5.78    | 0.61          | 1      |
| ellipf (Neg m)     | 0.87     | 0.80       | 2.79    | 3.75    | 0.54          | 1      |
| ellipeinc          | 1.04     | 0.81       | 9.01    | 24.00   | 2.99          | 1      |
| ellipeinc (Neg m)  | 1.24     | 0.99       | 4.25    | 5.76    | 1.08          | 1      |
| ellippiinc         | 52.78    | 0.78       | 1358.22 | 5821.91 | 1.52e5        | 50     |
| ellippiinc (Neg m) | 3.59     | 0.72       | 56.71   | 227.43  | 231.19        | 50     |
| ellippiinc (p.v.)  | 21.90    | 2.70       | 513.84  | 1487.45 | 10474.99      | 50     |

## Bulirsch's Elliptic Integrals
Bulirsh's elliptic integrals are not natively implemented in Wolfram Engine.
Nevertheless, some of the integrals can be converted to their Legendre's counterpart,
which are available on Wolfram Engine. However, for `cel` and `el2`, their values cannot
be mapped simply. Hence, the reference values are generated using the functions
submitted by Jan Mangaldan on [Wolfram Function Repository](https://resources.wolframcloud.com/FunctionRepository/).
As for `cel2`, it is mapped to `cel` with p=1.
| Function   | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|------------|----------|------------|---------|---------|---------------|--------|
| cel        | 0.62     | 0.51       | 3.55    | 37.90   | 3.58          | 1      |
| cel (p.v.) | 1.17     | 0.53       | 20.06   | 37.64   | 15.72         | 1      |
| cel1       | 0.50     | 0.00       | 8.68    | 8.68    | 1.61          | 1      |
| cel2       | 0.36     | 0.00       | 1.59    | 1.59    | 0.18          | 1      |
| el1        | 0.76     | 0.69       | 2.70    | 3.32    | 0.49          | 1      |
| el2        | 2.10     | 0.97       | 23.19   | 68.50   | 21.70         | 1      |
| el3        | 1119.33  | 0.68       | 9941.79 | 1.55e5  | 1.44e8        | 50     |
| el3 (p.v.) | 24.78    | 15.03      | 260.03  | 295.19  | 1407.09       | 50     |

## Carlson's Symmetric Elliptic Integrals

| Function       | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|----------------|----------|------------|---------|---------|---------------|--------|
| elliprf        | 5.45e6   | 0.60       | 5.43e7  | 5.90e7  | 1.89e14       | 1      |
| elliprg        | 5.64e5   | 0.82       | 8.33e6  | 9.34e6  | 2.72e12       | 50     |
| elliprj        | 4.13e5   | 1.17       | 1.08e7  | 2.35e7  | 4.33e12       | 50     |
| elliprj (p.v.) | 2.84e12  | 1.23       | 1.27e14 | 1.51e14 | 3.02e26       | 50     |
| elliprc        | 0.76     | 0.54       | 7.87    | 11.92   | 1.82          | 1      |
| elliprc (p.v.) | 0.65     | 0.64       | 2.03    | 2.64    | 0.31          | 1      |
| elliprd        | 9.65     | 2.30       | 105.66  | 1203.31 | 1924.29       | 50     |

## Carlson's Symmetric Elliptic Integrals with Boost Math
For functions that require higher precision, the accuracy of `ellip` on f64 is tested against
the answers from Boost Math at double precision.

| Function       | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|----------------|----------|------------|---------|---------|---------------|--------|
| elliprf        | 0.37     | 0.00       | 1.41    | 1.75    | 0.19          | 1      |
| elliprg        | 0.42     | 0.00       | 3.18    | 4.24    | 0.34          | 50     |
| elliprj        | 0.79     | 0.65       | 5.95    | 7.42    | 1.16          | 50     |
| elliprj (p.v.) | 5.52e11  | 0.65       | 2.54e12 | 1.58e14 | 5.81e25       | 50     |
