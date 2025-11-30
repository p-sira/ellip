# Contributing to Ellip

Thank you for your interest in contributing to Ellip. This document provides guidelines and information for contributors to help maintain the quality and accuracy of Ellip.

## Types of Contribution

In general, we accept:
- Documentation enhancements
- Typos and grammatical corrections.
- Bug fixes
- Technical insights
- Performance and accuracy improvements
- New code related to elliptic integral functions
- Unit tests

If you are unsure whether your contribution will fall under any of these categories, please open an issue and discuss.

## Licensing
All contribution will be licensed under Ellip's license, https://github.com/p-sira/ellip/blob/main/LICENSE. If you did not write the code yourself, ensure the existing license is compatible and include the license information with the contributed files, or obtain permission from the original author to relicense the contributed code.

## Getting Started

Before contributing, you should:

1. Read through the existing documentation at [docs.rs/ellip](https://docs.rs/ellip)
2. Check existing issues and pull requests to avoid duplication
3. For major changes, open an issue first to discuss your proposal

### Project Setup

```sh
# Clone the repository
git clone https://github.com/p-sira/ellip.git
cd ellip

# Build the project
cargo build

# Run tests
cargo test

# Run benchmarks
cargo bench
```

## Documentation Format

The documentation should be written in a concise, professional, technical manner. For functions, the documentation must follow this format:

- A **one-line description** of what the function does.
- **Mathematical definition.** For formulas, format into a Unicode-style mathematical notation using [Diagon](https://github.com/ArthurSonzogni/Diagon).
- **Parameters** of the function: include their valid domains.
- **Domain** of the function: list under what conditions the function may fail.
- **Graph** rendered using the implemented function (optional): see [ellip-plot-graph/](https://github.com/p-sira/ellip/tree/main/ellip-plot-graph).
- **Special Cases** (optional).
- **Notes** (optional).
- **Related Functions** (optional).
- **Examples**: standalone code showing basic usage of the function. Use `assert_close` for testing.
- **References**: literatures from reputable sources used in the implementation and documentation. The reference should be textbooks, online textbooks, peer-reviewed articles, and industry-standard implementations, such as Boost, GSL, or Mathematica.

Below is the template for documentation:

```rust
/// Computes [Function]
///
/// ```text
/// [Insert mathematical definition]
/// ```
///
/// ## Parameters
/// - 
///
/// ## Domain
/// - Returns error when:
///   - 
/// - Returns the Cauchy principal value if 
///
/// ## Graph
///
/// ## Special Cases
/// - 
///
/// ## Notes
///
/// # Related Functions
/// - 
///
/// # Examples
/// ```
/// use ellip::{[function], util::assert_close};
///
/// assert_close(...);
/// ```
///
/// # References
/// - 
```

For an example, see [ellip::ellippinc](https://docs.rs/ellip/latest/ellip/legendre/fn.ellippiinc.html).

## Code Contributions

### Style Guidelines

- Follow standard Rust format. Use `cargo fmt`.
- Run `cargo clippy` and address warnings.
- For mathematical variables, prefer standard notation.
- Add comments for non-obvious mathematical transformations.
- The documentation must follow the [guideline](#documentation-format).

### Function Signature

- Use functional style (e.g., `ellip_func(1.0, 3.0)`) rather than trait style (e.g., `(1.0).ellip_func(3.0)`).
- The types must be generic [`Float`](https://docs.rs/num-traits/0.2.19/num_traits/float/trait.Float.html) when applicable.
- Use `Result` types for fallible operations.

### Error Handling

- Use `Result` types for fallible operations.
- Error message must include the function name, followed by an informative description.
- Document conditions under which functions may fail.

### Testing

Ellip maintains high standards for numerical accuracy. All contributions must include appropriate tests, using reference values from reliable sources, such as textbooks, peer-reviewed articles, or industry-standard implementations.

### Unit Tests

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::assert_close;

    #[test]
    fn test_your_function() {
        let result = your_function(0.5).unwrap();
        assert_close(result, expected_value, 1e-15);
    }
}
```

You are encourage to use macros for tests to reduce boilerplates. Existing macros are `compare_test_data_boost` and `compare_test_data_wolfram`.

## Questions?

If you have questions about contributing:

- Open an issue for discussion
- Check existing issues and documentation
- Reach out to maintainers

Thank you for helping improve Ellip!
