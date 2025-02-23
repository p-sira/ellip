<h1 align="center">
    <a href="https://github.com/p-sira/ellip/">
        <img src="https://github.com/p-sira/ellip/blob/main/logo/ellip-logo.svg?raw=true" alt="ELLIP" width="300">
    </a>
</h1>

<p align="center">
    <a href="https://opensource.org/license/BSD-3-clause">
        <img src="https://img.shields.io/badge/License-BSD--3--Clause-brightgreen.svg" alt="License">
    </a>
    <a href="https://crates.io/crates/ellip">
        <img src="https://img.shields.io/crates/v/ellip" alt="Crate">
    </a>
    <a href="https://docs.rs/ellip">
        <img src="https://img.shields.io/badge/Docs-docs.rs-blue" alt="Documentation">
    </a>
</p>

<big><p align="center"> 
Elliptic integrals for Rust 
</p></big>

```shell
>> cargo add ellip
```

## Features
- Legendre's complete integrals
    - `ellipk`: Complete elliptic integral of the first kind.
    - `ellipe`: Complete elliptic integral of the second kind.
    - `ellippi`: Complete elliptic integral of the third kind.
    - `ellipd`: Complete elliptic integral of Legendre's type.
- Legendre's incomplete integrals
    - `ellipf`: Incomplete elliptic integral of the first kind.
    - `ellipeinc`: Incomplete elliptic integral of the second kind.
    - `ellippiinc`: Incomplete elliptic integral of the third kind.
    - `ellipdinc`: Incomplete elliptic integral of Legendre's type.
- Bulirsch's integrals
    - `cel`: General complete elliptic integral in Bulirsch's form.
    - `cel1`: Complete elliptic integral of the first kind in Bulirsch's form.
    - `cel2`: Complete elliptic integral of the second kind in Bulirsch's form.
    - `el1`: Incomplete elliptic integral of the first kind in Bulirsch's form.
    - `el2`: Incomplete elliptic integral of the second kind in Bulirsch's form.
    - `el3`: Incomplete elliptic integral of the third kind in Bulirsch's form.
- Carlson's symmetric integrals
    - `elliprf`: Symmetric elliptic integral of the first kind.
    - `elliprg`: Symmetric elliptic integral of the second kind.
    - `elliprj`: Symmetric elliptic integral of the third kind.
    - `elliprc`: Degenerate elliptic integral of RF.
    - `elliprd`: Degenerate elliptic integral of the third kind.

---

Learn more at [docs.rs](https://docs.rs/ellip).