# Ellip Rayon
Parallelized elliptic integral computation for Rust based on [Ellip](https://github.com/p-sira/ellip).

## Installation
```shell
cargo add ellip-rayon
```

## Machine-specific Threshold
Ellip Rayon employs parallelization if the length of the arguments exceesds certain thresholds. These thresholds depend on the core count, cache size, and architecture. For the most efficiency, these thresholds should be tuned on the target machine. 

1. Generating benchmark data
```shell
cargo bench
```
The process should take about 30-40 minutes to complete, and the file `benches/par_threshold.md` will be created. This reports the threshold for each function.

2. Using the generated threshold
```shell
cargo run --example generate_threshold_code
```
This script automatically replaces the thresholds in the source code.

3. Adding locally compiled library
From your working directory, run
```shell
cargo add --path path/to/your/ellip-rayon
```
