//! Genus 2 ramified model divisor arithmetic.
//!
//! This module provides explicit formulas for divisor addition and doubling
//! on ramified genus 2 hyperelliptic curves (one point at infinity).
//!
//! Three variants are provided optimized for different field characteristics:
//!
//! - [`arbitrary`] - Arbitrary characteristic (general case)
//! - [`not_char2`] - Not characteristic 2 (simplified formulas when h(x) = 0)
//! - [`char2`] - Characteristic 2 (uses XOR operations)

pub mod arbitrary;
pub mod char2;
pub mod not_char2;

#[cfg(test)]
mod blackbox_tests;
#[cfg(test)]
mod whitebox_tests;

// Re-export the commonly used types from arbitrary
pub use arbitrary::CurveConstants;
pub use arbitrary::DivisorCoords;
