# divisor-arithmetic

[![Crates.io](https://img.shields.io/crates/v/divisor-arithmetic.svg)](https://crates.io/crates/divisor-arithmetic)
[![Documentation](https://docs.rs/divisor-arithmetic/badge.svg)](https://docs.rs/divisor-arithmetic)
[![CI](https://github.com/salindne/divisor-arithmetic/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/salindne/divisor-arithmetic/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Explicit formulas for divisor arithmetic on genus 2 and 3 hyperelliptic curves.

This crate provides highly optimized implementations of divisor addition and doubling operations on the Jacobian of hyperelliptic curves, based on the research and Magma implementations from [salindne/divisorArithmetic](https://github.com/salindne/divisorArithmetic).

## Features

- **Genus 2 Ramified Model**: Curves with one point at infinity
  - `y² + h(x)y = f(x)` where `deg(f) = 5`, `deg(h) ≤ 2`
  - Three optimized variants:
    - `arbitrary` - Arbitrary characteristic (general case)
    - `not_char2` - Not characteristic 2 (simplified formulas when `h(x) = 0`)
    - `char2` - Characteristic 2 (XOR-based operations)

- **Genus 2 Split Model**: Curves with two points at infinity
  - `y² + h(x)y = f(x)` where `deg(f) = 6`, `deg(h) ≤ 3`
  - Balanced divisors carry an integer balance weight, reduced with respect to a
    positive (`Vpl`) or negative (`Vn = −Vpl − h`) basis — both bases provided
    (`add_neg`/`double_neg`, `add_pos`/`double_pos`)
  - Three variants: `not_char2`, `arbitrary` (any characteristic), `char2`
  - All formulas are cross-checked against a generic Cantor reference
    implementation (`generic::split`)

- **Field Implementations**:
  - `PrimeField<P>` - Prime fields F_p for small primes
  - `BinaryExtField<K>` - Binary extension fields GF(2^k) for k ≤ 24

- **Generic Polynomial Algorithms**: Reference implementations using Cantor's algorithm, NUCOMP, and NUDUPLE

## Usage

Add to your `Cargo.toml`:

```toml
[dependencies]
divisor-arithmetic = "0.1"
```

### Example: Divisor Addition on a Genus 2 Curve

```rust
use divisor_arithmetic::field::PrimeField;
use divisor_arithmetic::g2::ramified::not_char2::{add, double, CurveConstants, DivisorCoords};

type F31 = PrimeField<31>;

// Define curve: y² = x⁵ + 2x³ + 3x² + 4x + 5 over F_31
let curve = CurveConstants {
    f3: F31::new(2),
    f2: F31::new(3),
    f1: F31::new(4),
    f0: F31::new(5),
};

// Create degree-2 divisors
let d1 = DivisorCoords::deg2(F31::new(1), F31::new(2), F31::new(3), F31::new(4));
let d2 = DivisorCoords::deg2(F31::new(5), F31::new(6), F31::new(7), F31::new(8));

// Compute D1 + D2
let sum = add(&d1, &d2, &curve);

// Compute 2*D1
let doubled = double(&d1, &curve);
```

### Example: Characteristic 2 Fields

```rust
use divisor_arithmetic::field::BinaryExtField;
use divisor_arithmetic::g2::ramified::char2::{add, double, CurveConstants, DivisorCoords};

type GF256 = BinaryExtField<8>;  // GF(2^8)

// Define curve over GF(2^8)
let curve = CurveConstants {
    f2: GF256::new(0x1A),
    f1: GF256::new(0x2B),
    f0: GF256::new(0x3C),
    h2: GF256::new(0),
    h1: GF256::new(0x4D),
    h0: GF256::new(0x5E),
};

let d1 = DivisorCoords::deg2(
    GF256::new(0x11), GF256::new(0x22),
    GF256::new(0x33), GF256::new(0x44),
);

let doubled = double(&d1, &curve);
```

## Performance

Benchmark results on Apple M1 (single core):

Ramified model:

| Operation | Field | Time |
|-----------|-------|------|
| deg2 + deg2 (not_char2) | F_65521 | ~149 ns |
| deg2 + deg2 (arbitrary) | F_65521 | ~213 ns |
| deg2 + deg2 (char2) | GF(2^8) | ~185 ns |
| deg2 + deg2 (char2) | GF(2^16) | ~600 ns |
| 2*deg2 (not_char2) | F_65521 | ~194 ns |

Split model (degree-2 balanced divisors, negative basis):

| Operation | Field | Time |
|-----------|-------|------|
| add (not_char2) | F_65521 | ~180 ns |
| add (arbitrary) | F_65521 | ~190 ns |
| add (char2) | GF(2^8) | ~200 ns |
| add (char2) | GF(2^16) | ~625 ns |
| double (not_char2) | F_8191 | ~155 ns |
| double (char2) | GF(2^16) | ~680 ns |

Run benchmarks with:
```bash
cargo bench
```

## Testing

Run all tests:
```bash
cargo test --release
```

Or use the test script for detailed output:
```bash
./scripts/test.sh
```

## References

This implementation is based on:

- **Original Magma Implementation**: [salindne/divisorArithmetic](https://github.com/salindne/divisorArithmetic)
- **Thesis**: Sebastian Lindner, "Explicit Formulas for Hyperelliptic Curve Arithmetic", University of Calgary, 2020

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

