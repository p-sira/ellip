# Testing

This report denotes how to reproduce the tests and presents the accuracy of the ellip crate using **symmetric relative error**, defined as

![Symmetric relative error](https://github.com/p-sira/ellip/blob/main/examples/symmetric_error_on_white.svg?raw=true)

Errors are expressed in units of machine epsilon (ε). This report is generated using `cargo run --example generate_error_report`.

## Test Datasets

Ellip uses three sources of reference values:
1. **Wolfram Engine** for accuracy reporting and unit tests
2. **Boost.Math test data** in unit tests
3. **Original literature** in the unit tests for Bulirsh's integrals

The [Boost dataset](https://github.com/boostorg/math/tree/develop/test/) is included in the repository at [tests/data/boost/](https://github.com/p-sira/ellip/blob/main/tests/data/boost/).

The accuracy report shown here uses only the Wolfram dataset, as it provides uniform coverage across the domain of each function. The Wolfram reference values are generated using the scripts in [tests/wolfram/](https://github.com/p-sira/ellip/blob/main/tests/wolfram/). Each script samples inputs across the valid domain of the corresponding function, avoiding singularities by stopping short of the branch limit **μ**. This ensures meaningful reference values and prevents non-informative cases, such as large outputs that diverges toward infinity, where the observed error is dominated by floating-point limits.

To generate the Wolfram test data:

```sh
cd tests/wolfram
wolframscript --file [test-script].wls
```

The test datasets are not distributed with the crate by default. You may generate them locally using the scripts above or download the precomputed data from [tests/data/wolfram/](https://github.com/p-sira/ellip/tree/main/tests/data/wolfram/).

## f64 Results

{{ENV}}

### Legendre's Complete Elliptic Integrals

{{LEGENDRE_COMPLETE}}

### Legendre's Incomplete Elliptic Integrals

{{LEGENDRE_INCOMPLETE}}

### Bulirsch's Elliptic Integrals
Bulirsh's elliptic integrals are not natively implemented in Wolfram Engine. Nevertheless, some of the integrals can be converted to their Legendre's counterpart, which are available on Wolfram Engine. However, for `cel` and `el2`, their values cannot be mapped simply. Hence, the reference values are generated using the functions submitted by Jan Mangaldan on [Wolfram Function Repository](https://resources.wolframcloud.com/FunctionRepository/). As for `cel2`, it is mapped to `cel` with p=1.

{{BULIRSCH}}

### Carlson's Symmetric Elliptic Integrals

{{CARLSON}}

Note that `elliprj` is numerically unstable in the principal value domain when the symmetric parameter values are close (smaller than 1,000 epsilons).

### Miscellaneous Functions

{{MISC}}

## f32 Results

{{ENV_F32}}

### Legendre's Complete Elliptic Integrals

{{LEGENDRE_COMPLETE_F32}}

### Legendre's Incomplete Elliptic Integrals

{{LEGENDRE_INCOMPLETE_F32}}

### Bulirsch's Elliptic Integrals

{{BULIRSCH_F32}}

### Carlson's Symmetric Elliptic Integrals

{{CARLSON_F32}}

### Miscellaneous Functions

{{MISC_F32}}