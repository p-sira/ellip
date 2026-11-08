# Test Summary

See the [full error report](tests/README.md) for detailed accuracy results and test reproducibility information.

Benchmark on AMD Ryzen 5 4600H with Radeon Graphics @4.0 GHz running `x86_64-unknown-linux-gnu rustc 1.99.1` using ellip v1.2.4 with `libm` at `f64` precision (ε≈2.22e-16).

### Legendre's Elliptic Integrals
| Function            | Median Error (ε) | Max Error (ε) | Mean Performance |
|---------------------|------------------|---------------|------------------|
| ellipk              | 0.00             | 20.95         | 13.8 ns          |
| ellipe              | 0.00             | 1.99          | 12.5 ns          |
| ellipf              | 0.00             | 7.47          | 90.9 ns          |
| ellipeinc           | 0.00             | 24.66         | 126.1 ns         |
| ellippi             | 0.00             | 21.08         | 125.7 ns         |
| ellippiinc          | 0.00             | 395.31        | 191.4 ns         |
| ellippiinc_bulirsch | 0.00             | 395.31        | 157.5 ns         |
| ellipd              | 0.00             | 2.64          | 29.9 ns          |
| ellipdinc           | 0.00             | 8.38          | 85.6 ns          |
| ellipke             | 0.00             | 20.95         | 18.0 ns          |

### Bulirsch's Elliptic Integrals
| Function | Median Error (ε) | Max Error (ε) | Mean Performance |
|----------|------------------|---------------|------------------|
| cel      | 0.63             | 36.94         | 31.1 ns          |
| cel1     | 0.00             | 8.68          | 11.7 ns          |
| cel2     | 0.00             | 3.47          | 21.7 ns          |
| cel3     | 0.00             | 6.88          | 26.9 ns          |
| el1      | 0.00             | 1.70          | 37.2 ns          |
| el2      | 0.00             | 74.60         | 53.2 ns          |
| el3      | 0.00             | 18.65         | 81.1 ns          |

### Carlson's Symmetric Integrals
| Function | Median Error (ε) | Max Error (ε) | Mean Performance |
|----------|------------------|---------------|------------------|
| elliprf  | 0.00             | 1.57          | 52.1 ns          |
| elliprg  | 0.00             | 1.83          | 82.2 ns          |
| elliprj  | 0.58             | 168.11        | 154.8 ns         |
| elliprc  | 0.00             | 2.82          | 23.7 ns          |
| elliprd  | 0.00             | 2.68          | 71.1 ns          |

### Miscellaneous Functions
| Function      | Median Error (ε) | Max Error (ε) | Mean Performance |
|---------------|------------------|---------------|------------------|
| jacobi_zeta   | 0.00             | 8.66          | 135.6 ns         |
| heuman_lambda | 0.00             | 2.77          | 256.8 ns         |
