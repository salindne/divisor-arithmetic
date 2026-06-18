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

- **Batched group law** (ramified `not_char2`): `add_batch` / `double_batch`
  amortize the single field inversion across a whole batch of independent
  operations via Montgomery's trick (`field::batch_invert`) — the same strategy
  smalljac uses for generic-group order computations

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

### Field-operation counts

Wall-clock comparisons across implementations are confounded by field size,
operation definition, and language, so the explicit-formula literature compares
**field-operation counts** instead (field-size independent). Counts below are
for the generic-branch degree-2 `add`/`double` (negative basis), measured by
running the actual formulas over an instrumented field
(`cargo test --release g2::split::op_counts -- --nocapture`):

| Operation | M (mul) | S (sqr) | I (inv) | A (add/sub/dbl) |
|-----------|--------:|--------:|--------:|----------------:|
| not_char2 — add    | 26 | 2 | 1 | 37 |
| not_char2 — double | 31 | 3 | 1 | 38 |
| arbitrary — add    | 30 | 1 | 1 | 36 |
| arbitrary — double | 38 | 2 | 1 | 44 |
| char2 — add        | 27 | 1 | 1 | 34 |
| char2 — double     | 29 | 2 | 1 | 31 |

Each group operation uses exactly **one** field inversion (affine formulas).
These can be compared directly to the per-formula counts in Lange, Erickson–
Jacobson–Stein (real genus 2), and Costello–Lauter.

### Wall-clock comparison with smalljac (scalar and batched)

[smalljac](https://math.mit.edu/~drew/smalljac.html) (Andrew Sutherland) is a
highly optimized C library whose `hecurve_g2_compose` / `hecurve_g2_square`
implement the genus-2 **imaginary** (ramified, `deg f = 5`) group law — the same
model as this crate's `g2::ramified::not_char2`. Both are compared two ways:

- **scalar** — one field inversion per group operation (this crate's plain
  `add`/`double`; smalljac with `ctx = NULL`);
- **batched** — one field inversion shared across a batch of `N = 1024`
  independent operations via Montgomery's trick (this crate's
  [`add_batch`/`double_batch`] + [`field::batch_invert`]; smalljac's
  `hecurve_ctx_t` state machine + `ff_parallel_invert`). This is the throughput
  metric that matters for the generic-group order computations smalljac targets.

All numbers below were measured **on the same machine** (Apple Silicon, single
core; smalljac v4.1.3 + ff_poly v1.2.7 ported to arm64), on the same depressed
monic quintic (`f₄ = 0` — required, else smalljac falls back to generic Cantor),
in ns per group operation:

| field | op | this crate, scalar | smalljac, scalar | this crate, batched | smalljac, batched |
|-------|----|------------------:|-----------------:|-------------------:|------------------:|
| p = 65521 (16-bit) | add    | 222 | 100 | 99  | 47 |
| p = 65521 (16-bit) | double | 245 | 105 | 115 | 54 |
| 56-bit prime       | add    | 673 | 190 | 348 | 48 |
| 56-bit prime       | double | 691 | 210 | 384 | 54 |

Batching helps both implementations (it removes the inversion): this crate's
56-bit add drops 673 → 348 ns (1.9×), smalljac's 190 → 48 ns (4×). But smalljac
is faster in every cell, and the reason is the **field-arithmetic layer**, not
the genus-2 formulas (the [operation counts](#field-operation-counts) — ~26 M,
1 I — are essentially the same). Measured field ops at the 56-bit prime:

| field operation | this crate (`PrimeField`) | smalljac (ff_poly) |
|-----------------|--------------------------:|-------------------:|
| multiply        | 8.2 ns  | 4.5 ns |
| inversion       | 529 ns  | 157 ns |
| amortized inversion (batch of 1024) | ~25 ns | ~10 ns |

smalljac uses single-word **Montgomery** arithmetic (no division in the
multiply; a tuned binary-GCD inversion), whereas this crate's `PrimeField` uses
schoolbook `u128`-multiply + hardware modulo and an extended-Euclidean inversion
— ~1.8× slower to multiply and ~3.4× slower to invert. Since the scalar group op
is inversion-dominated (529 of 673 ns at 56-bit), that gap is what the batched
path removes, and it is what a Montgomery `PrimeField` would close. The takeaway:
**the explicit formulas and the batched group law are sound and competitive in
operation count; closing the wall-clock gap to smalljac is a field-arithmetic
optimization (Montgomery reduction), not a formula one.**

The harness (scalar + batched + raw field ops) and arm64 build notes are in
[`benches/smalljac-compare/`](benches/smalljac-compare/); run this crate's side
with `cargo bench -- not_char2`.

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

