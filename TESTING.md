# Test Summary

See the [full error report](tests/README.md) for detailed accuracy results and test reproducibility information.

Benchmark on AMD Ryzen 5 4600H with Radeon Graphics @4.0 GHz running `x86_64-unknown-linux-gnu rustc 1.98.1` using ellip v1.2.0 with `libm` at `f64` precision (ε≈2.22e-16).

### Legendre's Elliptic Integrals
| Function            | Median Error (ε) | Max Error (ε) | Mean Performance |
|---------------------|------------------|---------------|------------------|
| ellipk              | 0.00             | 20.95         | 14.0 ns          |
| ellipe              | 0.00             | 3.00          | 12.2 ns          |
| ellipf              | 0.00             | 7.47          | 105.8 ns         |
| ellipeinc           | 0.00             | 24.66         | 145.7 ns         |
| ellippi             | 0.00             | 21.08         | 134.4 ns         |
| ellippiinc          | 0.00             | 395.31        | 206.1 ns         |
| ellippiinc_bulirsch | 0.00             | 395.31        | 162.8 ns         |
| ellipd              | 0.00             | 2.64          | 32.5 ns          |
| ellipdinc           | 0.00             | 8.38          | 106.6 ns         |
| ellipke             | 0.00             | 20.95         | 16.8 ns          |

### Bulirsch's Elliptic Integrals
| Function | Median Error (ε) | Max Error (ε) | Mean Performance |
|----------|------------------|---------------|------------------|
| cel      | 0.63             | 36.94         | 31.0 ns          |
| cel1     | 0.00             | 8.68          | 11.9 ns          |
| cel2     | 0.00             | 3.47          | 21.8 ns          |
| cel3     | 0.00             | 6.88          | 26.4 ns          |
| el1      | 0.00             | 1.70          | 36.4 ns          |
| el2      | 0.00             | 74.60         | 52.4 ns          |
| el3      | 0.00             | 18.65         | 82.6 ns          |

### Carlson's Symmetric Integrals
| Function | Median Error (ε) | Max Error (ε) | Mean Performance |
|----------|------------------|---------------|------------------|
| elliprf  | 0.00             | 1.57          | 51.9 ns          |
| elliprg  | 0.00             | 5.25          | 82.8 ns          |
| elliprj  | 0.57             | 136.97        | 156.5 ns         |
| elliprc  | 0.00             | 2.82          | 23.8 ns          |
| elliprd  | 0.00             | 3.56          | 75.6 ns          |

### Miscellaneous Functions
| Function      | Median Error (ε) | Max Error (ε) | Mean Performance |
|---------------|------------------|---------------|------------------|
| jacobi_zeta   | 0.00             | 8.66          | 140.8 ns         |
| heuman_lambda | 0.00             | 2.77          | 270.8 ns         |
