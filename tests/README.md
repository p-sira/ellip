# Testing
This report presents the accuracy of the ellip crate using [**symmetric relative error**](https://www.boost.org/doc/libs/1_88_0/libs/math/doc/html/math_toolkit/relative_error.html) metric. Errors are expressed in units of machine epsilon (ε). The test data spans the domain of each function up to **μ** to avoid approaching the function's limit. The reference values are computed using [**Wolfram Engine**](https://www.wolfram.com/engine/). You can find the scripts in the directory [tests/wolfram/](https://github.com/p-sira/ellip/blob/main/tests/wolfram/). 
This report is generated on x86_64-unknown-linux-gnu rustc 1.90.0 using ellip v0.5.3 at `f64` precision (ε=2.2204460492503131e-16).

## Legendre's Complete Elliptic Integrals

| Function        | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|-----------------|----------|------------|---------|---------|---------------|--------|
| ellipk          | 0.27     | 0.00       | 1.51    | 20.95   | 0.81          | 1      |
| ellipk (Neg m)  | 0.58     | 0.00       | 1.55    | 108.14  | 15.91         | 1      |
| ellipe          | 0.30     | 0.00       | 1.82    | 3.00    | 0.18          | 1      |
| ellipe (Neg m)  | 0.42     | 0.51       | 1.64    | 1.95    | 0.20          | 1      |
| ellippi         | 0.41     | 0.00       | 1.65    | 20.95   | 0.60          | 1      |
| ellippi (Neg m) | 0.43     | 0.51       | 1.82    | 2.96    | 0.22          | 1      |
| ellippi (p.v.)  | 0.77     | 0.70       | 2.60    | 36.35   | 2.79          | 50     |
| ellipd          | 0.57     | 0.58       | 1.96    | 2.43    | 0.26          | 1      |
| ellipd (Neg m)  | 0.59     | 0.61       | 1.93    | 2.64    | 0.27          | 1      |

## Legendre's Incomplete Elliptic Integrals

| Function                    | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|-----------------------------|----------|------------|---------|---------|---------------|--------|
| ellipf                      | 0.65     | 0.66       | 2.54    | 7.47    | 0.44          | 1      |
| ellipf (Neg m)              | 0.60     | 0.66       | 2.02    | 3.00    | 0.29          | 1      |
| ellipeinc                   | 0.93     | 0.69       | 7.79    | 24.66   | 2.67          | 1      |
| ellipeinc (Neg m)           | 0.83     | 0.73       | 2.81    | 3.96    | 0.53          | 1      |
| ellippiinc                  | 1.49     | 0.72       | 18.93   | 165.83  | 62.92         | 1      |
| ellippiinc (Neg m)          | 1.11     | 0.68       | 13.27   | 25.98   | 4.99          | 1      |
| ellippiinc (p.v.)           | 11.39    | 2.92       | 114.13  | 1.04e3  | 3.38e3        | 1      |
| ellippiinc_bulirsch         | 1.67     | 0.77       | 18.93   | 165.83  | 62.85         | 1      |
| ellippiinc_bulirsch (Neg m) | 1.07     | 0.80       | 6.01    | 15.94   | 2.46          | 1      |
| ellipdinc                   | 1.31     | 1.08       | 4.24    | 8.38    | 1.09          | 1      |
| ellipdinc (Neg m)           | 1.16     | 0.93       | 3.55    | 4.20    | 0.82          | 1      |

## Bulirsch's Elliptic Integrals
Bulirsh's elliptic integrals are not natively implemented in Wolfram Engine. Nevertheless, some of the integrals can be converted to their Legendre's counterpart, which are available on Wolfram Engine. However, for `cel` and `el2`, their values cannot be mapped simply. Hence, the reference values are generated using the functions submitted by Jan Mangaldan on [Wolfram Function Repository](https://resources.wolframcloud.com/FunctionRepository/). As for `cel2`, it is mapped to `cel` with p=1.

| Function   | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|------------|----------|------------|---------|---------|---------------|--------|
| cel        | 0.97     | 0.75       | 5.64    | 28.60   | 2.86          | 1      |
| cel (p.v.) | 1.37     | 0.74       | 19.25   | 36.94   | 14.93         | 1      |
| cel1       | 0.54     | 0.00       | 7.82    | 8.68    | 1.56          | 1      |
| cel2       | 0.70     | 0.62       | 3.37    | 3.97    | 0.54          | 1      |
| el1        | 0.37     | 0.00       | 1.09    | 1.70    | 0.16          | 1      |
| el2        | 1.69     | 0.71       | 20.69   | 74.60   | 25.26         | 1      |
| el3        | 1.78     | 0.66       | 19.82   | 53.21   | 23.26         | 50     |
| el3 (p.v.) | 1.56     | 0.81       | 13.81   | 16.54   | 7.99          | 50     |

## Carlson's Symmetric Elliptic Integrals

| Function               | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|------------------------|----------|------------|---------|---------|---------------|--------|
| elliprf                | 0.40     | 0.51       | 1.41    | 1.57    | 0.18          | 1      |
| elliprg                | 0.45     | 0.51       | 2.60    | 5.25    | 0.37          | 1      |
| elliprj                | 0.81     | 0.66       | 6.20    | 7.42    | 1.22          | 50     |
| elliprj (p.v.)         | 1.51e3   | 0.63       | 141.81  | 7.12e5  | 8.29e8        | 1e10   |
| elliprj (p.v., small*) | 0.00     | 0.00       | 0.00    | 0.00    | 0.00          | 50     |
| elliprc                | 0.31     | 0.00       | 1.20    | 1.96    | 0.15          | 1      |
| elliprc (p.v.)         | 0.50     | 0.56       | 1.89    | 2.82    | 0.24          | 1      |
| elliprd                | 0.60     | 0.63       | 2.24    | 6.25    | 0.36          | 50     |

*small: Use small argument values, close to the function's limit. Results compared with Boost Math implementation without promoting double, i.e, computed purely using `f64`.

Current implementation of `elliprj` is less numerically stable in p.v. cases, as seen by large errors in the non-small test cases. That said, Ellip's results are consistent with Boost Math when limited to same precision (See [tests/data/boost/carlson.cpp](https://github.com/p-sira/ellip/blob/main/tests/data/boost/carlson.cpp)). Since the function is convergent, such errors can be mitigated when Rust's `f128` is released.

## Miscellaneous Functions

| Function            | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|---------------------|----------|------------|---------|---------|---------------|--------|
| jacobi_zeta         | 1.81     | 1.21       | 7.31    | 9.14    | 3.16          | 1      |
| jacobi_zeta (Neg m) | 1.89     | 1.42       | 7.67    | 8.89    | 3.14          | 1      |
| heuman_lambda       | 0.47     | 0.54       | 1.83    | 2.86    | 0.25          | 1      |

## f32 Implementation

Generated on x86_64-unknown-linux-gnu rustc 1.90.0 using ellip v0.5.3 at `f32` precision (ε≈1.19e-7).

### Legendre's Complete Elliptic Integrals

| Function        | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|-----------------|----------|------------|---------|---------|---------------|--------|
| ellipk          | 0.26     | 0.00       | 1.80    | 12.46   | 0.49          | 1      |
| ellipk (Neg m)  | 0.34     | 0.00       | 1.29    | 1.72    | 0.17          | 1      |
| ellipe          | 0.20     | 0.00       | 1.00    | 1.92    | 0.14          | 1      |
| ellipe (Neg m)  | 0.35     | 0.00       | 1.48    | 1.61    | 0.17          | 1      |
| ellippi         | 0.52     | 0.53       | 2.46    | 4.22    | 0.37          | 1      |
| ellippi (Neg m) | 0.48     | 0.55       | 2.05    | 2.92    | 0.27          | 1      |
| ellippi (p.v.)  | 0.80     | 0.70       | 2.45    | 34.47   | 2.59          | 50     |
| ellipd          | 0.64     | 0.59       | 2.56    | 15.09   | 0.85          | 1      |
| ellipd (Neg m)  | 0.55     | 0.59       | 1.87    | 2.16    | 0.25          | 1      |

### Legendre's Incomplete Elliptic Integrals

| Function                    | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|-----------------------------|----------|------------|---------|---------|---------------|--------|
| ellipf                      | 0.47     | 0.54       | 1.97    | 5.20    | 0.29          | 1      |
| ellipf (Neg m)              | 0.46     | 0.53       | 1.60    | 3.54    | 0.23          | 1      |
| ellipeinc                   | 0.70     | 0.57       | 5.01    | 24.02   | 1.93          | 1      |
| ellipeinc (Neg m)           | 0.72     | 0.69       | 2.39    | 3.12    | 0.39          | 1      |
| ellippiinc                  | 1.05     | 0.68       | 15.59   | 71.79   | 11.40         | 1      |
| ellippiinc (Neg m)          | 1.15     | 0.67       | 12.36   | 40.74   | 8.56          | 1      |
| ellippiinc (p.v.)           | 5.24     | 1.65       | 55.19   | 298.48  | 327.85        | 1      |
| ellippiinc_bulirsch         | 1.47e6   | 0.79       | 1.68e7  | 1.68e7  | 2.25e13       | 1      |
| ellippiinc_bulirsch (Neg m) | 2.25e6   | 0.85       | 1.68e7  | 1.68e7  | 3.26e13       | 1      |
| ellipdinc                   | 0.92     | 0.74       | 3.75    | 11.97   | 1.00          | 1      |
| ellipdinc (Neg m)           | 0.84     | 0.76       | 2.92    | 4.44    | 0.52          | 1      |

### Bulirsch's Elliptic Integrals

| Function   | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|------------|----------|------------|---------|---------|---------------|--------|
| cel        | 0.85     | 0.67       | 10.30   | 18.59   | 2.40          | 1      |
| cel (p.v.) | 0.98     | 0.73       | 9.17    | 14.06   | 2.50          | 1      |
| cel1       | 0.37     | 0.00       | 1.58    | 1.59    | 0.18          | 1      |
| cel2       | 1.27e5   | 0.64       | 8.39e6  | 8.39e6  | 1.05e12       | 1      |
| el1        | 0.41     | 0.00       | 1.93    | 1.96    | 0.23          | 1      |
| el2        | 1.51     | 0.69       | 20.14   | 91.18   | 19.39         | 1      |
| el3        | 1.09     | 0.65       | 8.17    | 17.86   | 2.52          | 50     |
| el3 (p.v.) | 1.40     | 0.86       | 22.73   | 23.50   | 7.49          | 50     |

### Carlson's Symmetric Elliptic Integrals

| Function               | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|------------------------|----------|------------|---------|---------|---------------|--------|
| elliprf                | 0.38     | 0.51       | 1.63    | 1.79    | 0.18          | 1      |
| elliprg                | 0.36     | 0.00       | 1.45    | 2.12    | 0.18          | 1      |
| elliprj                | 0.80     | 0.62       | 7.74    | 14.15   | 1.89          | 50     |
| elliprj (p.v.)         | 1.94     | 0.70       | 34.47   | 116.46  | 64.04         | 1e10   |
| elliprj (p.v., small*) | NAN      | NAN        | NAN     | NAN     | NAN           | 50     |
| elliprc                | 0.30     | 0.00       | 1.14    | 1.73    | 0.14          | 1      |
| elliprc (p.v.)         | 0.52     | 0.57       | 1.81    | 2.15    | 0.23          | 1      |
| elliprd                | 0.52     | 0.58       | 1.85    | 5.82    | 0.30          | 50     |

### Miscellaneous Functions

| Function            | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|---------------------|----------|------------|---------|---------|---------------|--------|
| jacobi_zeta         | 1.40     | 0.97       | 6.16    | 7.32    | 1.89          | 1      |
| jacobi_zeta (Neg m) | 1.52     | 1.17       | 5.75    | 6.71    | 1.82          | 1      |
| heuman_lambda       | 0.51     | 0.56       | 2.01    | 2.80    | 0.29          | 1      |