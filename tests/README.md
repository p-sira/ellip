# Testing
This report presents the accuracy of the ellip crate using [**symmetric relative error**](https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/relative_error.html) metric. Errors are expressed in units of machine epsilon (ε). The test data spans the domain of each function up to **μ** to avoid approaching the function's limit. The reference values are computed using [**Wolfram Engine**](https://www.wolfram.com/engine/). You can find the scripts in the directory [tests/wolfram/](https://github.com/p-sira/ellip/blob/main/tests/wolfram/). 
This report is generated on x86_64-unknown-linux-gnu rustc 1.89.0 using ellip v0.3.7 at `f64` precision (ε≈2.22e-16).

## Legendre's Complete Elliptic Integrals

| Function        | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|-----------------|----------|------------|---------|---------|---------------|--------|
| ellipk          | 0.27     | 0.00       | 1.52    | 20.95   | 0.81          | 1      |
| ellipk (Neg m)  | 0.58     | 0.00       | 1.71    | 108.14  | 15.91         | 1      |
| ellipe          | 0.30     | 0.00       | 1.83    | 3.00    | 0.18          | 1      |
| ellipe (Neg m)  | 0.42     | 0.51       | 1.71    | 1.95    | 0.20          | 1      |
| ellippi         | 0.41     | 0.00       | 1.68    | 20.95   | 0.60          | 1      |
| ellippi (Neg m) | 0.43     | 0.51       | 1.83    | 2.96    | 0.22          | 1      |
| ellippi (p.v.)  | 0.77     | 0.70       | 2.66    | 36.35   | 2.79          | 50     |
| ellipd          | 0.57     | 0.58       | 2.03    | 2.43    | 0.26          | 1      |
| ellipd (Neg m)  | 0.59     | 0.61       | 1.93    | 2.64    | 0.27          | 1      |

## Legendre's Incomplete Elliptic Integrals

| Function           | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|--------------------|----------|------------|---------|---------|---------------|--------|
| ellipf             | 0.65     | 0.66       | 2.59    | 7.47    | 0.44          | 1      |
| ellipf (Neg m)     | 0.60     | 0.66       | 2.04    | 3.00    | 0.29          | 1      |
| ellipeinc          | 0.93     | 0.69       | 8.15    | 24.66   | 2.67          | 1      |
| ellipeinc (Neg m)  | 0.83     | 0.73       | 2.90    | 3.96    | 0.53          | 1      |
| ellippiinc         | 1.49     | 0.72       | 22.57   | 165.83  | 62.92         | 1      |
| ellippiinc (Neg m) | 1.11     | 0.68       | 14.53   | 25.98   | 4.99          | 1      |
| ellippiinc (p.v.)  | 11.39    | 2.92       | 168.79  | 1.04e3  | 3.38e3        | 1      |
| ellipdinc          | 1.31     | 1.08       | 4.30    | 8.38    | 1.09          | 1      |
| ellipdinc (Neg m)  | 1.16     | 0.93       | 3.55    | 4.20    | 0.82          | 1      |

## Bulirsch's Elliptic Integrals
Bulirsh's elliptic integrals are not natively implemented in Wolfram Engine. Nevertheless, some of the integrals can be converted to their Legendre's counterpart, which are available on Wolfram Engine. However, for `cel` and `el2`, their values cannot be mapped simply. Hence, the reference values are generated using the functions submitted by Jan Mangaldan on [Wolfram Function Repository](https://resources.wolframcloud.com/FunctionRepository/). As for `cel2`, it is mapped to `cel` with p=1.

| Function   | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|------------|----------|------------|---------|---------|---------------|--------|
| cel        | 0.92     | 0.71       | 5.64    | 28.60   | 2.83          | 1      |
| cel (p.v.) | 1.33     | 0.70       | 20.06   | 38.34   | 15.92         | 1      |
| cel1       | 0.54     | 0.00       | 8.68    | 8.68    | 1.56          | 1      |
| cel2       | 0.61     | 0.61       | 2.27    | 3.97    | 0.40          | 1      |
| el1        | 0.36     | 0.00       | 1.11    | 1.60    | 0.15          | 1      |
| el2        | 1.60     | 0.70       | 18.01   | 79.92   | 24.32         | 1      |
| el3        | 1.66     | 0.66       | 19.40   | 46.32   | 15.89         | 50     |
| el3 (p.v.) | 1.56     | 0.81       | 13.88   | 16.54   | 7.99          | 50     |

## Carlson's Symmetric Elliptic Integrals

| Function               | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|------------------------|----------|------------|---------|---------|---------------|--------|
| elliprf                | 0.37     | 0.00       | 1.41    | 1.75    | 0.19          | 1      |
| elliprg                | 0.39     | 0.00       | 2.25    | 2.45    | 0.25          | 1      |
| elliprj                | 0.80     | 0.67       | 6.59    | 7.42    | 1.15          | 50     |
| elliprj (p.v.)         | 2.04e3   | 0.69       | 87.56   | 7.60e5  | 1.07e9        | 1e10   |
| elliprj (p.v., small*) | 0.00     | 0.00       | 0.00    | 0.00    | 0.00          | 50     |
| elliprc                | 0.30     | 0.00       | 1.23    | 1.96    | 0.15          | 1      |
| elliprc (p.v.)         | 0.49     | 0.54       | 1.91    | 2.82    | 0.25          | 1      |
| elliprd                | 0.60     | 0.62       | 2.41    | 6.49    | 0.37          | 50     |

*small: Use small argument values, close to the function's limit. Results compared with Boost Math implementation without promoting double, i.e, computed purely using `f64`.

Current implementation of `elliprj` is less numerically stable in p.v. cases, as seen by large errors in the non-small test cases. That said, Ellip's results are consistent with Boost Math when limited to same precision (See [tests/data/boost/carlson.cpp](https://github.com/p-sira/ellip/blob/main/tests/data/boost/carlson.cpp)). Since the function is convergent, such errors can be mitigated when Rust's `f128` is released.
