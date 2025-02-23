# Changelog
## 0.3
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
- `cel`: Bulirsch's general complete elliptic integral
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