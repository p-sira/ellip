# Changelog
## 0.3
### 0.3.0
**New Features**
- All functions now support generic Float.
- `unchecked`: module for direct access to functions without argument checking and special case evaluation
- `ellippi`: Complete elliptic integral of the third kind.

**Bug Fixes**
- `elliprf`: incorrect result for RF(x,x,0)

**Changes**
- Add dependency `num-traits`.

## 0.2
### 0.2.1
**Improvements**
- Reduce crate size by removing the logo and test data from the source. To perform full test, download the test data from the [repository](https://github.com/p-sira/ellip/tree/main/tests/data).

### 0.2.0
**New Features**
- `cel`: Bulirsch's General complete elliptic integral
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
- Fix domain error message in `elliprd`

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
- Add CHANGELOG

### 0.1.0
**New Features**
- `ellipf`: Incomplete elliptic integral of the first kind.
- `ellipeinc`: Incomplete elliptic integral of the second kind.
- `ellipk`: Complete elliptic integral of the first kind.
- `ellipe`: Complete elliptic integral of the second kind.
- `elliprf`: Symmetric elliptic integral of the first kind.
- `elliprd`: Degenerate elliptic integral of the third kind.
- `elliprc`: Degenerate elliptic integral of RF.