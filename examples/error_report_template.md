# Testing
This report presents the accuracy of the ellip crate using **symmetric relative error**, defined as

![Symmetric relative error](https://github.com/p-sira/ellip/blob/main/examples/symmetric_error.svg?raw=true)

Errors are expressed in units of machine epsilon (ε). The test data spans the domain of each function up to **μ** to avoid approaching the function's limit. The reference values are computed using [**Wolfram Engine**](https://www.wolfram.com/engine/). You can find the scripts in the directory [tests/wolfram/](https://github.com/p-sira/ellip/blob/main/tests/wolfram/). 
{{ENV}}

## Legendre's Complete Elliptic Integrals

{{LEGENDRE_COMPLETE}}

## Legendre's Incomplete Elliptic Integrals

{{LEGENDRE_INCOMPLETE}}

## Bulirsch's Elliptic Integrals
Bulirsh's elliptic integrals are not natively implemented in Wolfram Engine. Nevertheless, some of the integrals can be converted to their Legendre's counterpart, which are available on Wolfram Engine. However, for `cel` and `el2`, their values cannot be mapped simply. Hence, the reference values are generated using the functions submitted by Jan Mangaldan on [Wolfram Function Repository](https://resources.wolframcloud.com/FunctionRepository/). As for `cel2`, it is mapped to `cel` with p=1.

{{BULIRSCH}}

## Carlson's Symmetric Elliptic Integrals

{{CARLSON}}

Note that `elliprj` is numerically unstable in the principal value domain when the symmetric parameter values are close (smaller than 1,000 epsilons).

## Miscellaneous Functions

{{MISC}}

## f32 Implementation

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