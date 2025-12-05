#!/bin/bash
set -e

# =============================================================================
# Divisor Arithmetic Test Suite
# =============================================================================
# This script runs all whitebox and blackbox tests for genus 2 ramified formulas.
# Exit on first failure.

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
YELLOW='\033[1;33m'
BOLD='\033[1m'
NC='\033[0m'

print_header() {
    echo ""
    echo -e "${BLUE}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║${NC} ${BOLD}$1${NC}"
    echo -e "${BLUE}╚══════════════════════════════════════════════════════════════╝${NC}"
}

print_section() {
    echo -e "\n${CYAN}▶${NC} ${BOLD}$1${NC}"
}

print_success() {
    echo -e "${GREEN}✓${NC} $1"
}

print_fail() {
    echo -e "${RED}✗${NC} $1"
}

# Change to project root (parent of scripts/)
cd "$(dirname "$0")/.."

print_header "Divisor Arithmetic Test Suite"
echo -e "${YELLOW}Testing genus 2 ramified model formulas${NC}"
echo ""

# Track results
TOTAL_PASSED=0
TOTAL_FAILED=0

# Helper to extract test count from output
extract_count() {
    echo "$1" | grep "test result:" | sed 's/.*ok\. \([0-9]*\) passed.*/\1/' | head -1
}

# -----------------------------------------------------------------------------
# Field Tests
# -----------------------------------------------------------------------------
print_section "Field Arithmetic Tests (PrimeField, BinaryExtField)"
OUTPUT=$(cargo test field::tests --release 2>&1) || true
if echo "$OUTPUT" | grep -q "test result: ok"; then
    FIELD_COUNT=$(extract_count "$OUTPUT")
    print_success "Field tests passed (${FIELD_COUNT:-?} tests)"
    TOTAL_PASSED=$((TOTAL_PASSED + ${FIELD_COUNT:-0}))
else
    print_fail "Field tests failed"
    TOTAL_FAILED=$((TOTAL_FAILED + 1))
fi

# -----------------------------------------------------------------------------
# Whitebox Tests
# -----------------------------------------------------------------------------
print_section "Whitebox Tests (66 cases: 22 arbitrary + 22 char2 + 22 not_char2)"

# ARBITRARY whitebox
OUTPUT=$(cargo test g2::ramified::whitebox_tests::arbitrary_tests --release 2>&1) || true
if echo "$OUTPUT" | grep -q "test result: ok"; then
    ARB_COUNT=$(extract_count "$OUTPUT")
    print_success "ARBITRARY whitebox tests passed (${ARB_COUNT:-22} tests)"
    TOTAL_PASSED=$((TOTAL_PASSED + ${ARB_COUNT:-22}))
else
    print_fail "ARBITRARY whitebox tests failed"
    TOTAL_FAILED=$((TOTAL_FAILED + 1))
fi

# CHAR2 whitebox
OUTPUT=$(cargo test g2::ramified::whitebox_tests::char2_tests --release 2>&1) || true
if echo "$OUTPUT" | grep -q "test result: ok"; then
    CH2_COUNT=$(extract_count "$OUTPUT")
    print_success "CHAR2 whitebox tests passed (${CH2_COUNT:-22} tests)"
    TOTAL_PASSED=$((TOTAL_PASSED + ${CH2_COUNT:-22}))
else
    print_fail "CHAR2 whitebox tests failed"
    TOTAL_FAILED=$((TOTAL_FAILED + 1))
fi

# NOT_CHAR2 whitebox
OUTPUT=$(cargo test g2::ramified::whitebox_tests::not_char2_tests --release 2>&1) || true
if echo "$OUTPUT" | grep -q "test result: ok"; then
    NCH2_COUNT=$(extract_count "$OUTPUT")
    print_success "NOT_CHAR2 whitebox tests passed (${NCH2_COUNT:-22} tests)"
    TOTAL_PASSED=$((TOTAL_PASSED + ${NCH2_COUNT:-22}))
else
    print_fail "NOT_CHAR2 whitebox tests failed"
    TOTAL_FAILED=$((TOTAL_FAILED + 1))
fi

# -----------------------------------------------------------------------------
# Blackbox Tests
# -----------------------------------------------------------------------------
print_section "Blackbox Tests (11 tests across field sizes)"

# ARBITRARY blackbox
OUTPUT=$(cargo test g2::ramified::blackbox_tests::arbitrary_blackbox --release 2>&1) || true
if echo "$OUTPUT" | grep -q "test result: ok"; then
    ARB_BB_COUNT=$(extract_count "$OUTPUT")
    print_success "ARBITRARY blackbox tests passed (${ARB_BB_COUNT:-4} tests)"
    TOTAL_PASSED=$((TOTAL_PASSED + ${ARB_BB_COUNT:-4}))
else
    print_fail "ARBITRARY blackbox tests failed"
    TOTAL_FAILED=$((TOTAL_FAILED + 1))
fi

# CHAR2 blackbox
OUTPUT=$(cargo test g2::ramified::blackbox_tests::char2_blackbox --release 2>&1) || true
if echo "$OUTPUT" | grep -q "test result: ok"; then
    CH2_BB_COUNT=$(extract_count "$OUTPUT")
    print_success "CHAR2 blackbox tests passed (${CH2_BB_COUNT:-3} tests)"
    TOTAL_PASSED=$((TOTAL_PASSED + ${CH2_BB_COUNT:-3}))
else
    print_fail "CHAR2 blackbox tests failed"
    TOTAL_FAILED=$((TOTAL_FAILED + 1))
fi

# NOT_CHAR2 blackbox
OUTPUT=$(cargo test g2::ramified::blackbox_tests::not_char2_blackbox --release 2>&1) || true
if echo "$OUTPUT" | grep -q "test result: ok"; then
    NCH2_BB_COUNT=$(extract_count "$OUTPUT")
    print_success "NOT_CHAR2 blackbox tests passed (${NCH2_BB_COUNT:-4} tests)"
    TOTAL_PASSED=$((TOTAL_PASSED + ${NCH2_BB_COUNT:-4}))
else
    print_fail "NOT_CHAR2 blackbox tests failed"
    TOTAL_FAILED=$((TOTAL_FAILED + 1))
fi

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
if [ $TOTAL_FAILED -eq 0 ]; then
    echo -e "${GREEN}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${GREEN}║                    All tests passed!                         ║${NC}"
    printf "${GREEN}║           Total: %-4d tests passed                           ║${NC}\n" $TOTAL_PASSED
    echo -e "${GREEN}╚══════════════════════════════════════════════════════════════╝${NC}"
else
    echo -e "${RED}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${RED}║                  Some tests failed!                          ║${NC}"
    printf "${RED}║         Passed: %-4d  Failed: %-4d                          ║${NC}\n" $TOTAL_PASSED $TOTAL_FAILED
    echo -e "${RED}╚══════════════════════════════════════════════════════════════╝${NC}"
    exit 1
fi
echo ""
