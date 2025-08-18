# Changelog
## 0.3
### 0.3.7
**Bug Fixes**
- `ellippiinc`: Fix negative m domain where n < 0 and m <= n.

**Improvements**
- Improve performance across all functions by removing unnecessary checks.
- Improve and fix error messages.
- Add test summary and benchmark data in README.
- `el3`: Improve accuracy.
- `elliprj`: Slightly increase accuracy of the p.v. cases.

**Tests**
- Add Wolfram test data for `ellipd` and `ellipdinc`.

### 0.3.6
**Bug Fixes**
- `ellipdinc`: Fix incorrect result when m = -inf.
- `el3`: Fix NAN when p is small and kc = 0.

**Improvements**
- Improve performance of `legendre` functions by lazy checking and `el3` by using array instead of Vec.
- Handle special cases and edge cases (infinity, nan, etc.) for all `carlson` functions.
- Increase test coverage.

### 0.3.5
**Deprecations**
- Deprecate `BulirschConst`: The trait signature will be changed from `BulirschConst` to `BulirschConst<T>` and the trait will be relocated under `unstable` flag.

**Bug Fixes**
- `cel1`: Fix incorrect values in negative kc cases.
- `el3`: Fix incorrect values in `pm > 0.5` branch.

**Improvements**
- Loosen the bounds of the `bulirsch` functions.
- Improve the performance of the `bulirsch` functions by lazy checking.
- Handle special cases and edge cases (infinity, nan, etc.) for all `bulirsch` functions.
- Add special cases documentation for `bulirsch` functions.
- Improve and fix documentation.
- Reduce package size.

**Tests**
- Increase test coverage.
- Include `bulirsch` functions in the benchmark module.
- Use native Wolfram functions to generate test data for `cel1`, `el1`, and `el3`.
- Add error report for `cel1` and `cel2`.  

### 0.3.4
**Bug Fixes**
- `ellippiinc`: Fix not returning error when n sin²φ = 1.

**Improvements**
- Handle special cases and edge cases (infinity, nan, etc.) for all `legendre` functions.
- Add special cases documentation for `legendre` functions.
- Increase test coverage.
- Report code coverage.
- Move examples into separate crates, reducing dev dependencies and test time.

### 0.3.3
**Changes**
- Use `numeric_literals` to parse floats within the functions to make the code more readable.

### 0.3.2
**Improvements**
- Reduce package size.
- Use custom error type `StrErr` instead of `&'static str` to improve readability.
- Format documentation.

### 0.3.1
**Improvements**
- Improve documentation and add graphs and examples for each function.
- Improve the precision of the `ellipk`, `ellipf`, `ellipe`, `ellipeinc`, `ellippiinc`, and `elliprj` functions.
- Improve the speed of `ellippiinc` function.
- `el3`: Add error when parameter goes out of the function's range (|kc| > 0 for p < 0) and improved the error message when 1 + px² equals zero.

**New Features**
- Add `assert_close` function in `util` module.

**Bug Fixes**
- `ellipd`: Fix incorrect answers for m < 0. The function also returns infinity instead of throwing error at m = 1.
- Fix the domains for `ellipf` and `ellipeinc`.
- `el3`: Fix index out of bound in small el3 cases.
- `ellippiinc`: Fix infinite loop leading to stack overflow.

### 0.3.0
**New Features**
- `ellippi`: Complete elliptic integral of the third kind.
- `ellippiinc`: Incomplete elliptic integral of the third kind.
- `ellipd`: Complete elliptic integral of Legendre's type.
- `ellipdinc`: Incomplete elliptic integral of Legendre's type.
- `cel1`: Complete elliptic integral of the first kind in Bulirsch's form.
- `cel2`: Complete elliptic integral of the second kind in Bulirsch's form.
- `el1`: Incomplete elliptic integral of the first kind in Bulirsch's form.
- `el2`: Incomplete elliptic integral of the second kind in Bulirsch's form.
- `el3`: Incomplete elliptic integral of the third kind in Bulirsch's form.
- `BulirschConst` trait controls the precision of cel and el functions.
- All functions now support generic Float.

**Bug Fixes**
- `elliprf`: Fix incorrect result for RF(x,x,0).
- `elliprj` and `elliprc`: Fix losing precision due to ln(1+x).
- `ellipe`: Fix incorrect result during internal condition |t| > 10.

**Changes**
- Add dependency `num-traits` and `num-lazy`.

## 0.2
### 0.2.1
**Improvements**
- Reduce crate size by removing the logo and test data from the source. To perform full test, download the test data from the [repository](https://github.com/p-sira/ellip/tree/main/tests/data).

### 0.2.0
**New Features**
- `cel`: Bulirsch's general complete elliptic integral.
- `elliprg`: Symmetric elliptic integral of the second kind.
- `elliprj`: Symmetric elliptic integral of the third kind.

**Improvements**
- `elliprf`: Add special cases handling to improve performance and accuracy.
- `elliprd`: Use Boost Math implementation to improve accuracy.
- `ellipe`: Explicitly make internal functions inline.

**Changes**
- Change error messages in `ellipk`, `ellipf`, `ellipe`, and `ellipeinc`, to describe without mathematical notations.
- Split error message in `elliprd`.

**Bug Fixes**
- Fix domain error message in `elliprd`.

**Others**
- Add Ellip logo
- Add links to references in the documentations. 

## 0.1
### 0.1.2
**Bug Fixes**
- Fix `ellipe` logic in negative m cases.

### 0.1.1
**Others**
- Update README to reflect the crate's functionalities.
- Add CHANGELOG.

### 0.1.0
**New Features**
- `ellipf`: Incomplete elliptic integral of the first kind.
- `ellipeinc`: Incomplete elliptic integral of the second kind.
- `ellipk`: Complete elliptic integral of the first kind.
- `ellipe`: Complete elliptic integral of the second kind.
- `elliprf`: Symmetric elliptic integral of the first kind.
- `elliprd`: Degenerate elliptic integral of the third kind.
- `elliprc`: Degenerate elliptic integral of RF.