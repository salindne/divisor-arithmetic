//! Genus 2 **split** model divisor arithmetic (two points at infinity).
//!
//! Split genus 2 curves have the form `y² + h(x)·y = f(x)` with `deg(f) = 6`,
//! `deg(h) ≤ 3`. Divisors use the balanced 4-coordinate Mumford representation
//! `(u, v, w, n)` with an integer balance weight `n`, reduced with respect to a
//! positive (`Vpl`) or negative (`Vn = −Vpl − h`) basis.
//!
//! Variants (matching the ramified module's layout, one per characteristic):
//! - [`not_char2`] — fields of characteristic ≠ 2 (`h = 0`)
//!
//! Each variant provides `precompute` (curve constants + reduced basis) plus
//! `add_neg`/`double_neg` and `add_pos`/`double_pos`.

pub mod not_char2;

/// Record that a named formula branch (matching the Magma `DEBUG` labels, e.g.
/// `"DBL03"`) was taken. Compiles to a no-op outside tests.
#[cfg(not(test))]
#[inline(always)]
pub(crate) fn branch(_label: &'static str) {}

/// Test build: record taken branches so tests can assert full coverage.
#[cfg(test)]
#[inline]
pub(crate) fn branch(label: &'static str) {
    coverage::hit(label);
}

#[cfg(test)]
pub(crate) mod coverage {
    use std::collections::BTreeSet;
    use std::sync::{Mutex, OnceLock};

    fn store() -> &'static Mutex<BTreeSet<&'static str>> {
        static HITS: OnceLock<Mutex<BTreeSet<&'static str>>> = OnceLock::new();
        HITS.get_or_init(|| Mutex::new(BTreeSet::new()))
    }

    pub fn hit(label: &'static str) {
        store().lock().unwrap().insert(label);
    }

    /// Labels (from `expected`) that have not been hit yet.
    pub fn missing<'a>(expected: &[&'a str]) -> Vec<&'a str> {
        let hits = store().lock().unwrap();
        expected.iter().copied().filter(|l| !hits.contains(*l)).collect()
    }
}

#[cfg(test)]
mod test_support;
#[cfg(test)]
mod wb_vectors;
#[cfg(test)]
mod whitebox_tests;
#[cfg(test)]
mod blackbox_tests;
