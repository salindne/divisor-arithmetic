//! Genus 2 divisor arithmetic.
//!
//! This module provides explicit formulas for divisor addition and doubling
//! on genus 2 hyperelliptic curves.
//!
//! ## Curve Models
//!
//! - **Ramified Model**: Curves with one point at infinity
//!   - `y² + h(x)y = f(x)` where `deg(f) = 5`, `deg(h) ≤ 2`
//!   - See [`ramified`] module with variants for different characteristics:
//!     - [`ramified::arbitrary`] - Arbitrary characteristic (general case)
//!     - [`ramified::not_char2`] - Not characteristic 2 (simplified, h = 0)
//!     - [`ramified::char2`] - Characteristic 2
//!
//! - **Split Model**: Curves with two points at infinity  
//!   - `y² + h(x)y = f(x)` where `deg(f) = 6`, `deg(h) ≤ 3`

pub mod ramified;
// pub mod split;  // TODO: implement split model

#[cfg(test)]
mod tests; // Integration tests comparing specialized vs generic
