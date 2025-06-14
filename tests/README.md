# Testing
This report presents the accuracy of the ellip crate using [**symmetric relative error**](https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/relative_error.html)
metric. Errors are expressed in units of machine epsilon (ε).
The test data spans the domain of each function up to **μ** to avoid approaching the function's limit.
The reference values are computed using [**Wolfram Engine**](https://www.wolfram.com/engine/).
You can find the scripts in the directory [tests/wolfram/](https://github.com/p-sira/ellip/blob/main/tests/wolfram/).
This report is generated on x86_64-unknown-linux-gnu rustc 1.87.0 using ellip v0.3.1 at `f64` precision (ε≈2.22e-16).

## Legendre's Complete Elliptic Integrals
| Function        | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|-----------------|----------|------------|---------|---------|---------------|--------|
| ellipk          | 0.51     | 0.00       | 1.88    | 131.90  | 18.15         | 1      |
| ellipk (Neg m)  | 0.70     | 0.52       | 2.20    | 107.50  | 15.78         | 1      |
| ellipe          | 0.55     | 0.67       | 2.29    | 3.00    | 0.30          | 1      |
| ellipe (Neg m)  | 0.66     | 0.75       | 2.24    | 3.03    | 0.35          | 1      |
| ellippi         | 24.39    | 0.00       | 881.35  | 1216.37 | 16048.11      | 50     |
| ellippi (Neg m) | 65.04    | 0.55       | 2470.49 | 3197.93 | 129979.02     | 50     |
| ellippi (p.v.)  | 3.36     | 0.71       | 45.99   | 520.23  | 582.48        | 50     |

## Legendre's Incomplete Elliptic Integrals
| Function           | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|--------------------|----------|------------|---------|---------|---------------|--------|
| ellipf             | 0.81     | 0.72       | 3.27    | 5.78    | 0.61          | 1      |
| ellipf (Neg m)     | 0.87     | 0.80       | 2.79    | 3.75    | 0.54          | 1      |
| ellipeinc          | 1.04     | 0.81       | 9.01    | 24.00   | 2.99          | 1      |
| ellipeinc (Neg m)  | 1.24     | 0.99       | 4.25    | 5.76    | 1.08          | 1      |
| ellippiinc         | 52.78    | 0.78       | 1358.22 | 5821.91 | 151961.95     | 50     |
| ellippiinc (Neg m) | 3.59     | 0.72       | 56.71   | 227.43  | 231.19        | 50     |
| ellippiinc (p.v.)  | 21.90    | 2.70       | 513.84  | 1487.45 | 10474.99      | 50     |

## Bulirsch's Elliptic Integrals
For Bulirsh's elliptic integrals, the reference values are generated using the function
submitted by Jan Mangaldan on [Wolfram Function Repository](https://resources.wolframcloud.com/FunctionRepository/).
| Function   | Mean (ε) | Median (ε) | P99 (ε) | Max (ε)  | Variance (ε²) | μ (ε²) |
|------------|----------|------------|---------|----------|---------------|--------|
| cel        | 0.62     | 0.51       | 3.55    | 37.90    | 3.58          | 1      |
| cel (p.v.) | 1.17     | 0.53       | 20.06   | 37.64    | 15.72         | 1      |
| el1        | 1.01     | 0.82       | 3.92    | 12.22    | 1.24          | 1      |
| el2        | 2.10     | 0.97       | 23.19   | 68.50    | 21.70         | 1      |
| el3        | 27.18    | 2.61       | 143.89  | 11585.44 | 223942.40     | 50     |
| el3 (p.v.) | 31.67    | 20.26      | 260.08  | 295.19   | 2060.17       | 50     |

## Bulirsch's Incomplete Elliptic Integrals

## Carlson's Symmetric Elliptic Integrals

