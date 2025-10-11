# Ellip Rayon
Parallelized elliptic integral computation for Rust based on [Ellip](https://github.com/p-sira/ellip).

## Machine-specific Threshold
Ellip Rayon employs parallelization if the length of the arguments reaches the thresholds. For the most efficiency, these thresholds should be tuned on the target machine. By running `cargo bench`, the tuning should take about 30-40 minutes to complete and the file `benches/par_threshold.md` will be created. This reports the threshold for each function, which you can replace them in the sourcecode.