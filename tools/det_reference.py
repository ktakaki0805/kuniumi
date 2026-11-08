#!/usr/bin/env python3
"""Correctly rounded reference values for `kuniumi_math::det` (ticket M1-01).

Writes the Rust module `crates/kuniumi-math/src/det/golden.rs` to standard
output. Run from the repository root (see tools/README.md):

    tools/.venv/bin/python tools/det_reference.py \
        > crates/kuniumi-math/src/det/golden.rs

Method:

- Inputs are fixed in advance, below. Python floats are IEEE 754 binary64 and
  parse decimal literals with correct rounding, like Rust.
- Each value is computed with mpmath at 256 and at 512 bits of precision. Both
  must round to the same double, otherwise the script fails.
- The 512-bit value is converted to an exact `Fraction` and rounded to the
  nearest double (ties to even, subnormals included) by Python's exact
  `Fraction` -> `float` conversion. This is the correctly rounded result.
- The script also records on which side of that double the exact value lies
  (+1 above, -1 below, 0 if the result is exact). The two doubles around the
  exact value are the only faithfully rounded results (error below 1 ULP).
  libm (FreeBSD msun) is not correctly rounded, so at some inputs the tests
  expect the other one of the two; see `det/tests.rs`.
"""

from __future__ import annotations

import math
import struct
import sys
from fractions import Fraction

import mpmath
from mpmath import mp, mpf

PRECISIONS_BITS = (256, 512)

SIN_COS_INPUTS = [1e-5, 0.5, -0.7, 1.0, 2.0, 3.0, math.pi, 10.0, 100.0, 1e5, 1e22]
TAN_INPUTS = [1e-5, 0.5, 1.0, -1.2, 1.5, 2.0, math.pi / 2, 10.0, 1e3, 1e5, 1e22]
ATAN_INPUTS = [1e-8, 0.1, 0.5, -0.75, 1.0, 1.5, 3.0, 10.0, 1e3, -1e8]
ASIN_INPUTS = [-1e-8, 0.1, 0.25, 0.3, 0.5, -0.5, 0.7, 0.9, 0.99, 0.999999]
ACOS_INPUTS = [1e-8, 0.1, 0.25, 0.5, -0.5, 0.7, 0.9, 0.99, -0.99, 0.999999]
EXP_INPUTS = [1e-10, 0.5, 1.0, -1.0, 2.0, 10.0, -10.0, 100.0, 700.0, -700.0, -740.0]
LN_INPUTS = [
    5e-324,
    1e-10,
    0.1,
    0.5,
    0.9999999,
    1.0000001,
    2.0,
    3.0,
    10.0,
    123456.789,
    1e10,
    1e300,
]

# Inputs for the two-argument functions: atan2(y, x) and pow(x, y).
ATAN2_INPUTS = [
    (1.0, 1.0),
    (1.0, -1.0),
    (-1.0, -1.0),
    (-1.0, 1.0),
    (0.5, 2.0),
    (2.0, 0.5),
    (3.0, -4.0),
    (-1e-10, 1.0),
    (1e10, 1e-10),
    (1.0, -1e10),
]
POW_INPUTS = [
    (2.0, 0.5),
    (0.5, 0.5),
    (2.5, 3.7),
    (123.456, 2.5),
    (7.0, 1.0 / 3.0),
    (0.3, -50.0),
    (1e6, 50.0),
    (10.0, -3.0),
    (1.0000001, 1e7),
    (3.0, 0.1),
]

FUNCTIONS_1 = [
    ("SIN", "sin", mpmath.sin, SIN_COS_INPUTS),
    ("COS", "cos", mpmath.cos, SIN_COS_INPUTS),
    ("TAN", "tan", mpmath.tan, TAN_INPUTS),
    ("ATAN", "atan", mpmath.atan, ATAN_INPUTS),
    ("ASIN", "asin", mpmath.asin, ASIN_INPUTS),
    ("ACOS", "acos", mpmath.acos, ACOS_INPUTS),
    ("EXP", "exp", mpmath.exp, EXP_INPUTS),
    ("LN", "ln", mpmath.log, LN_INPUTS),
]
FUNCTIONS_2 = [
    ("ATAN2", "atan2", "y, x", mpmath.atan2, ATAN2_INPUTS),
    ("POW", "pow", "x, y", mpmath.power, POW_INPUTS),
]
CONSTANTS = [
    ("PI", "π", lambda: mpmath.pi),
    ("FRAC_PI_2", "π/2", lambda: mpmath.pi / 2),
    ("FRAC_PI_4", "π/4", lambda: mpmath.pi / 4),
    ("FRAC_3PI_4", "3π/4", lambda: 3 * mpmath.pi / 4),
]


def bits(x: float) -> int:
    """The IEEE 754 binary64 bit pattern of `x`."""
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def rust_hex(b: int) -> str:
    """A u64 as a Rust hex literal with 4-digit groups."""
    h = f"{b:016X}"
    return "0x" + "_".join(h[i : i + 4] for i in range(0, 16, 4))


def to_fraction(v: mpf) -> Fraction:
    """The exact rational value of a finite mpf."""
    if not mpmath.isfinite(v):
        raise ValueError(f"non-finite reference value {v}")
    numerator, denominator = v.as_integer_ratio()
    return Fraction(int(numerator), int(denominator))


def round_to_double(exact: Fraction) -> tuple[float, int, Fraction]:
    """Nearest double to `exact`, the side of the exact value, and the distance.

    The side is +1 if the exact value is above the double, -1 below, 0 equal.
    The distance is in ULPs, measured with the spacing between the double and
    its neighbour on the side of the exact value (which matters at powers of
    two).
    """
    r = float(exact)  # exact conversion with round-to-nearest-even
    if Fraction(r) == exact:
        return r, 0, Fraction(0)
    side = 1 if exact > Fraction(r) else -1
    neighbour = math.nextafter(r, math.inf if side > 0 else -math.inf)
    spacing = abs(Fraction(neighbour) - Fraction(r))
    return r, side, abs(exact - Fraction(r)) / spacing


def reference(fn, args: tuple[float, ...]) -> tuple[float, int, Fraction]:
    """Correctly rounded fn(*args), checked at two precisions."""
    results = []
    for prec in PRECISIONS_BITS:
        with mp.workprec(prec):
            value = fn(*(mpf(a) for a in args))
            results.append(round_to_double(to_fraction(value)))
    (r_lo, side_lo, _), (r_hi, side_hi, d_hi) = results
    if bits(r_lo) != bits(r_hi) or side_lo != side_hi:
        raise SystemExit(f"precision disagreement for {args}: {r_lo!r} vs {r_hi!r}")
    return r_hi, side_hi, d_hi


def self_check() -> None:
    """Checks the conversion helpers against values Python already knows."""
    with mp.workprec(PRECISIONS_BITS[-1]):
        for x in [0.5, -0.7, 1e22, -1e-300, 5e-324, -5e-324, 1.7976931348623157e308]:
            r, side, d = round_to_double(to_fraction(mpf(x)))
            assert bits(r) == bits(x) and side == 0 and d == 0, x
        # math.pi and math.e are the correctly rounded values of π and e.
        r, side, _ = round_to_double(to_fraction(mpmath.pi))
        assert bits(r) == bits(math.pi) and side == 1  # π = 3.14159265358979323846…
        r, side, _ = round_to_double(to_fraction(-mpmath.pi / 2))
        assert bits(r) == bits(-math.pi / 2) and side == -1
        r, side, _ = round_to_double(to_fraction(mpmath.e))
        assert bits(r) == bits(math.e) and side == 1  # e = 2.71828182845904523536…


def side_word(side: int) -> str:
    return {1: "above", -1: "below", 0: "equal to"}[side]


def main() -> None:
    self_check()
    out = []
    w = out.append
    w("//! Correctly rounded reference values for the `det` golden tests.")
    w("//!")
    w(f"//! Generated by `tools/det_reference.py` with mpmath {mpmath.__version__};")
    w("//! do not edit by hand. Each function's entries are")
    w("//! `(input bits…, correctly rounded result bits, side)`, where `side` is")
    w("//! +1 if the exact result lies above the rounded one, −1 below, 0 if the")
    w("//! result is exact. The comments give the distance from the exact value.")
    for name, label, thunk in CONSTANTS:
        w("")
        with mp.workprec(PRECISIONS_BITS[-1]):
            r, _, d = round_to_double(to_fraction(thunk()))
        w(f"/// {label} rounded to the nearest double ({float(d):.3f} ULP away).")
        w(f"pub(super) const {name}: u64 = {rust_hex(bits(r))};")
    for const, fname, fn, inputs in FUNCTIONS_1:
        w("")
        w(f"/// `{fname}`: (x, {fname}(x) correctly rounded, side of the exact value).")
        w(f"pub(super) const {const}: &[(u64, u64, i8)] = &[")
        for x in inputs:
            r, side, d = reference(fn, (x,))
            w(f"    // {fname}({x!r}) = {r!r}, exact value {float(d):.3f} ULP {side_word(side)}")
            w(f"    ({rust_hex(bits(x))}, {rust_hex(bits(r))}, {side}),")
        w("];")
    for const, fname, argnames, fn, inputs in FUNCTIONS_2:
        w("")
        w(
            f"/// `{fname}`: ({argnames}, {fname}({argnames}) correctly rounded, "
            "side of the exact value)."
        )
        w(f"pub(super) const {const}: &[(u64, u64, u64, i8)] = &[")
        for a, b in inputs:
            r, side, d = reference(fn, (a, b))
            w(
                f"    // {fname}({a!r}, {b!r}) = {r!r}, "
                f"exact value {float(d):.3f} ULP {side_word(side)}"
            )
            # One value per line, the layout rustfmt gives these tuples.
            w("    (")
            for v in (a, b, r):
                w(f"        {rust_hex(bits(v))},")
            w(f"        {side},")
            w("    ),")
        w("];")
    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
