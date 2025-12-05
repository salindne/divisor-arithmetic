//! # Divisor Arithmetic
//!
//! Explicit formulas for divisor arithmetic on hyperelliptic curves.
//!
//! This library implements efficient polynomial-level formulas for adding and
//! doubling divisors on hyperelliptic curves of genus 2 and 3.
//!
//! ## Curve Models
//!
//! - **Ramified Model**: Curves with one point at infinity (genus 2)
//!   - `y² + h(x)y = f(x)` where `deg(f) = 5`, `deg(h) ≤ 2`
//!
//! - **Split Model**: Curves with two points at infinity (genus 2 and 3)
//!   - `y² + h(x)y = f(x)` where `deg(f) = 2g+2`, `deg(h) ≤ g+1`
//!
//! ## Usage
//!
//! ```rust,ignore
//! use divisor_arithmetic::field::PrimeField;
//! use divisor_arithmetic::g2::ramified::{add, double};
//! ```

pub mod field;
pub mod poly;

pub mod g2;
pub mod g3;
pub mod generic;

/// Coefficient extraction helper - mirrors Magma's Coeff function
#[inline]
pub fn coeff<T: Clone + Default>(coeffs: &[T], n: usize) -> T {
    if n < coeffs.len() {
        coeffs[n].clone()
    } else {
        T::default()
    }
}
