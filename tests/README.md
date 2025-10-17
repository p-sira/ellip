# Testing
This report presents the accuracy of the ellip crate using **symmetric relative error**, defined as

![Symmetric relative error](https://github.com/p-sira/ellip/blob/main/examples/symmetric_error.svg?raw=true)

Errors are expressed in units of machine epsilon (ε). The test data spans the domain of each function up to **μ** to avoid approaching the function's limit. The reference values are computed using [**Wolfram Engine**](https://www.wolfram.com/engine/). You can find the scripts in the directory [tests/wolfram/](https://github.com/p-sira/ellip/blob/main/tests/wolfram/). 
This report is generated on x86_64-unknown-linux-gnu rustc 1.90.0 using ellip v0.5.6 at `f64` precision (ε=2.2204460492503131e-16).

## Legendre's Complete Elliptic Integrals

| Function        | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|-----------------|----------|------------|---------|---------|---------------|--------|
| ellipk          | 0.27     | 0.00       | 1.51    | 20.95   | 0.81          | 1      |
| ellipk (Neg m)  | 0.54     | 0.00       | 1.55    | 108.14  | 15.93         | 1      |
| ellipe          | 0.30     | 0.00       | 1.82    | 3.00    | 0.18          | 1      |
| ellipe (Neg m)  | 0.42     | 0.51       | 1.64    | 1.95    | 0.20          | 1      |
| ellippi         | 0.41     | 0.00       | 1.65    | 20.95   | 0.60          | 1      |
| ellippi (Neg m) | 0.43     | 0.51       | 1.82    | 2.96    | 0.22          | 1      |
| ellippi (p.v.)  | 0.17     | 0.00       | 2.11    | 36.35   | 2.56          | 50     |
| ellipd          | 0.44     | 0.00       | 1.96    | 2.43    | 0.32          | 1      |
| ellipd (Neg m)  | 0.10     | 0.00       | 1.87    | 2.64    | 0.16          | 1      |

## Legendre's Incomplete Elliptic Integrals

| Function                    | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|-----------------------------|----------|------------|---------|---------|---------------|--------|
| ellipf                      | 0.33     | 0.00       | 2.39    | 7.47    | 0.42          | 1      |
| ellipf (Neg m)              | 0.28     | 0.00       | 1.95    | 3.00    | 0.26          | 1      |
| ellipeinc                   | 0.51     | 0.00       | 7.50    | 24.66   | 2.69          | 1      |
| ellipeinc (Neg m)           | 0.50     | 0.00       | 2.68    | 3.38    | 0.54          | 1      |
| ellippiinc                  | 0.89     | 0.00       | 14.58   | 165.83  | 49.26         | 1      |
| ellippiinc (Neg m)          | 0.59     | 0.00       | 8.53    | 25.98   | 3.87          | 1      |
| ellippiinc (p.v.)           | 7.27     | 1.76       | 110.88  | 395.31  | 738.07        | 1      |
| ellippiinc_bulirsch         | 1.01     | 0.00       | 14.58   | 165.83  | 49.39         | 1      |
| ellippiinc_bulirsch (Neg m) | 0.50     | 0.00       | 5.14    | 15.94   | 1.69          | 1      |
| ellipdinc                   | 0.21     | 0.00       | 3.62    | 8.38    | 0.56          | 1      |
| ellipdinc (Neg m)           | 0.12     | 0.00       | 2.98    | 4.20    | 0.29          | 1      |

## Bulirsch's Elliptic Integrals
Bulirsh's elliptic integrals are not natively implemented in Wolfram Engine. Nevertheless, some of the integrals can be converted to their Legendre's counterpart, which are available on Wolfram Engine. However, for `cel` and `el2`, their values cannot be mapped simply. Hence, the reference values are generated using the functions submitted by Jan Mangaldan on [Wolfram Function Repository](https://resources.wolframcloud.com/FunctionRepository/). As for `cel2`, it is mapped to `cel` with p=1.

| Function   | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|------------|----------|------------|---------|---------|---------------|--------|
| cel        | 0.75     | 0.65       | 3.47    | 28.60   | 2.04          | 1      |
| cel (p.v.) | 1.09     | 0.59       | 19.25   | 36.94   | 13.66         | 1      |
| cel1       | 0.54     | 0.00       | 7.82    | 8.68    | 1.56          | 1      |
| cel2       | 0.50     | 0.00       | 2.73    | 3.47    | 0.51          | 1      |
| el1        | 0.09     | 0.00       | 1.09    | 1.70    | 0.08          | 1      |
| el2        | 0.37     | 0.00       | 5.54    | 74.60   | 16.74         | 1      |
| el3        | 1.15     | 0.00       | 19.40   | 53.21   | 20.34         | 50     |
| el3 (p.v.) | 1.24     | 0.00       | 13.81   | 16.54   | 8.60          | 50     |

## Carlson's Symmetric Elliptic Integrals

| Function       | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|----------------|----------|------------|---------|---------|---------------|--------|
| elliprf        | 0.33     | 0.00       | 1.41    | 1.57    | 0.19          | 1      |
| elliprg        | 0.26     | 0.00       | 2.60    | 5.25    | 0.38          | 1      |
| elliprj        | 0.75     | 0.58       | 6.20    | 7.42    | 1.27          | 50     |
| elliprj (p.v.) | 1.07     | 0.51       | 12.05   | 136.97  | 28.82         | 1000   |
| elliprc        | 0.20     | 0.00       | 1.20    | 1.96    | 0.14          | 1      |
| elliprc (p.v.) | 0.14     | 0.00       | 1.75    | 2.82    | 0.15          | 1      |
| elliprd        | 0.49     | 0.00       | 2.24    | 6.25    | 0.40          | 50     |

Note that `elliprj` is numerically unstable in the principal value domain when the symmetric parameter values are close (smaller than 1,000 epsilons).

## Miscellaneous Functions

| Function            | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|---------------------|----------|------------|---------|---------|---------------|--------|
| jacobi_zeta         | 0.07     | 0.00       | 3.12    | 7.59    | 0.32          | 1      |
| jacobi_zeta (Neg m) | 0.07     | 0.00       | 3.87    | 8.66    | 0.36          | 1      |
| heuman_lambda       | 0.37     | 0.00       | 1.82    | 2.86    | 0.24          | 1      |

## f32 Implementation

Generated on x86_64-unknown-linux-gnu rustc 1.90.0 using ellip v0.5.6 at `f32` precision (ε≈1.19e-7).

### Legendre's Complete Elliptic Integrals

| Function        | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|-----------------|----------|------------|---------|---------|---------------|--------|
| ellipk          | 0.26     | 0.00       | 1.80    | 12.46   | 0.49          | 1      |
| ellipk (Neg m)  | 0.30     | 0.00       | 1.29    | 1.72    | 0.18          | 1      |
| ellipe          | 0.20     | 0.00       | 1.00    | 1.92    | 0.14          | 1      |
| ellipe (Neg m)  | 0.35     | 0.00       | 1.48    | 1.61    | 0.17          | 1      |
| ellippi         | 0.52     | 0.53       | 2.46    | 4.22    | 0.37          | 1      |
| ellippi (Neg m) | 0.48     | 0.55       | 2.02    | 2.92    | 0.26          | 1      |
| ellippi (p.v.)  | 0.18     | 0.00       | 2.14    | 34.47   | 2.36          | 50     |
| ellipd          | 0.51     | 0.00       | 2.56    | 15.09   | 0.92          | 1      |
| ellipd (Neg m)  | 0.07     | 0.00       | 1.83    | 2.16    | 0.11          | 1      |

### Legendre's Incomplete Elliptic Integrals

| Function                    | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|-----------------------------|----------|------------|---------|---------|---------------|--------|
| ellipf                      | 0.23     | 0.00       | 1.95    | 5.20    | 0.25          | 1      |
| ellipf (Neg m)              | 0.20     | 0.00       | 1.56    | 3.54    | 0.19          | 1      |
| ellipeinc                   | 0.37     | 0.00       | 4.84    | 24.02   | 1.89          | 1      |
| ellipeinc (Neg m)           | 0.39     | 0.00       | 2.31    | 3.12    | 0.39          | 1      |
| ellippiinc                  | 0.53     | 0.00       | 8.72    | 71.78   | 7.85          | 1      |
| ellippiinc (Neg m)          | 0.64     | 0.00       | 11.93   | 40.74   | 6.10          | 1      |
| ellippiinc (p.v.)           | 2.97     | 0.00       | 42.49   | 128.05  | 99.10         | 1      |
| ellippiinc_bulirsch         | 1.16e6   | 0.00       | 1.68e7  | 1.68e7  | 1.82e13       | 1      |
| ellippiinc_bulirsch (Neg m) | 1.86e6   | 0.00       | 1.68e7  | 1.68e7  | 2.78e13       | 1      |
| ellipdinc                   | 0.21     | 0.00       | 3.59    | 11.97   | 0.69          | 1      |
| ellipdinc (Neg m)           | 0.10     | 0.00       | 2.13    | 4.44    | 0.21          | 1      |

### Bulirsch's Elliptic Integrals

| Function   | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|------------|----------|------------|---------|---------|---------------|--------|
| cel        | 0.55     | 0.00       | 2.44    | 9.41    | 0.59          | 1      |
| cel (p.v.) | 0.70     | 0.56       | 6.48    | 13.37   | 1.54          | 1      |
| cel1       | 0.37     | 0.00       | 1.58    | 1.59    | 0.18          | 1      |
| cel2       | 0.41     | 0.00       | 2.00    | 2.00    | 0.33          | 1      |
| el1        | 0.14     | 0.00       | 1.36    | 1.89    | 0.13          | 1      |
| el2        | 0.19     | 0.00       | 2.58    | 91.18   | 9.95          | 1      |
| el3        | 0.59     | 0.00       | 7.58    | 17.86   | 2.43          | 50     |
| el3 (p.v.) | 1.06     | 0.00       | 22.73   | 23.50   | 7.98          | 50     |

### Carlson's Symmetric Elliptic Integrals

| Function       | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|----------------|----------|------------|---------|---------|---------------|--------|
| elliprf        | 0.31     | 0.00       | 1.63    | 1.79    | 0.19          | 1      |
| elliprg        | 0.17     | 0.00       | 1.45    | 2.12    | 0.15          | 1      |
| elliprj        | 0.70     | 0.56       | 6.84    | 14.15   | 1.62          | 50     |
| elliprj (p.v.) | 115.48   | 0.53       | 11.34   | 1.54e5  | 1.46e7        | 1000   |
| elliprc        | 0.19     | 0.00       | 1.14    | 1.73    | 0.13          | 1      |
| elliprc (p.v.) | 0.17     | 0.00       | 1.83    | 2.15    | 0.17          | 1      |
| elliprd        | 0.41     | 0.00       | 1.85    | 5.82    | 0.32          | 50     |

### Miscellaneous Functions

| Function            | Mean (ε) | Median (ε) | P99 (ε) | Max (ε) | Variance (ε²) | μ (ε²) |
|---------------------|----------|------------|---------|---------|---------------|--------|
| jacobi_zeta         | 0.02     | 0.00       | 0.42    | 3.99    | 0.06          | 1      |
| jacobi_zeta (Neg m) | 0.02     | 0.00       | 0.00    | 4.52    | 0.08          | 1      |
| heuman_lambda       | 0.40     | 0.00       | 1.98    | 2.80    | 0.29          | 1      |