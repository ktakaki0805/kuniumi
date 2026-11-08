#!/usr/bin/env python3
"""Error of `kuniumi_math::det` against correctly rounded values (ticket M1-01).

Generates seeded random inputs for every function, evaluates them with the
example `crates/kuniumi-math/examples/det_eval.rs` (through cargo) and compares
each result with the exact value computed by mpmath. Prints, per function, the
maximum error in ULPs of the exact value, the share of correctly rounded
results and the worst input. The numbers are quoted in the documentation of
`kuniumi_math::det`. Run from the repository root (see tools/README.md):

    mise exec -- tools/.venv/bin/python tools/det_error.py [count]

Method:

- The error of a result y is |y - f(x)| / ulp(f(x)), where ulp(v) is the
  spacing of doubles in the binade of v (2^-1074 for subnormals).
- f(x) is computed with mpmath at 160 bits. The 20 largest errors of each
  function are recomputed at 600 bits and must agree to 1e-9 ULP.
- Input ranges follow the accuracy tests of `det` (det/tests.rs), plus y and x
  uniform in [-10, 10] for atan2.
"""

from __future__ import annotations

import math
import random
import struct
import subprocess
import sys
from fractions import Fraction

import mpmath
from mpmath import mp, mpf

SEED = 20261108
WORK_PREC_BITS = 160
CHECK_PREC_BITS = 600
CHECK_WORST = 20

MPMATH = {
    "sin": mpmath.sin,
    "cos": mpmath.cos,
    "tan": mpmath.tan,
    "atan": mpmath.atan,
    "asin": mpmath.asin,
    "acos": mpmath.acos,
    "exp": mpmath.exp,
    "ln": mpmath.log,
    "atan2": mpmath.atan2,
    "pow": mpmath.power,
}


def bits(x: float) -> int:
    return struct.unpack("<Q", struct.pack("<d", x))[0]


def from_bits(b: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", b))[0]


def ordered(x: float) -> int:
    b = bits(x)
    return -(b & 0x7FFF_FFFF_FFFF_FFFF) if b >> 63 else b


def from_ordered(i: int) -> float:
    return from_bits((1 << 63 | -i) if i < 0 else i)


def uniform(rng: random.Random, n: int, lo: float, hi: float) -> list[float]:
    return [lo + (hi - lo) * rng.random() for _ in range(n)]


def log_uniform(rng: random.Random, n: int, lo: float, hi: float) -> list[float]:
    """Positive doubles with uniform bit patterns between lo and hi."""
    return [from_bits(rng.randint(bits(lo), bits(hi))) for _ in range(n)]


def signed(rng: random.Random, xs: list[float]) -> list[float]:
    return [-x if rng.random() < 0.5 else x for x in xs]


def near_half_pi_multiples(rng: random.Random, n: int) -> list[float]:
    out = []
    for _ in range(n):
        k = rng.randint(-63_661, 63_661)
        out.append(from_ordered(ordered(k * (math.pi / 2)) + rng.randint(-(1 << 16), 1 << 16)))
    return out


def random_normal(rng: random.Random, biased_exponent: int) -> float:
    return from_bits(
        rng.getrandbits(1) << 63 | biased_exponent << 52 | rng.getrandbits(52)
    )


def inputs(rng: random.Random, n: int) -> dict[str, list[tuple[float, ...]]]:
    """Input tuples per function; n per function."""
    tenth = n // 10
    trig = lambda: uniform(rng, n - tenth, -1e5, 1e5) + near_half_pi_multiples(rng, tenth)
    third = n // 3
    atan2 = list(zip(uniform(rng, third, -10, 10), uniform(rng, third, -10, 10)))
    for i in range(n - third):
        ey = rng.randint(1, 2046)
        ex = rng.randint(1, 2046) if i % 2 == 0 else min(max(ey + rng.randint(-64, 64), 1), 2046)
        atan2.append((random_normal(rng, ey), random_normal(rng, ex)))
    half = n // 2
    bases = [1e6 * (1.0 - rng.random()) for _ in range(half)]
    bases += log_uniform(rng, n - half, 2.0**-20, 1e6)
    return {
        "sin": [(x,) for x in trig()],
        "cos": [(x,) for x in trig()],
        "tan": [(x,) for x in trig()],
        "atan": [
            (x,)
            for x in uniform(rng, half, -1e5, 1e5)
            + signed(rng, log_uniform(rng, n - half, 2.0**-60, 2.0**60))
        ],
        "asin": [(x,) for x in uniform(rng, n, -1.0, 1.0)],
        "acos": [(x,) for x in uniform(rng, n, -1.0, 1.0)],
        "exp": [(x,) for x in uniform(rng, n, -745.0, 709.0)],
        "ln": [
            (x,)
            for x in log_uniform(rng, half, 5e-324, sys.float_info.max)
            + uniform(rng, n - half, 0.5, 2.0)
        ],
        "atan2": atan2,
        "pow": [(x, rng.uniform(-50.0, 50.0)) for x in bases],
    }


def evaluate(cases: dict[str, list[tuple[float, ...]]]) -> dict[str, list[float]]:
    """Results of det for every input, from the det_eval example."""
    lines = []
    for name, args in cases.items():
        lines += [" ".join([name] + [f"{bits(a):x}" for a in arg]) for arg in args]
    run = subprocess.run(
        ["cargo", "run", "--release", "--quiet", "--locked", "-p", "kuniumi-math"]
        + ["--example", "det_eval"],
        input="\n".join(lines) + "\n",
        capture_output=True,
        text=True,
        check=True,
    )
    values = iter(from_bits(int(w, 16)) for w in run.stdout.split())
    return {name: [next(values) for _ in args] for name, args in cases.items()}


def ulp_error(y: float, exact: mpf) -> Fraction:
    """|y - exact| in ULPs of the binade of exact."""
    ex = Fraction(*(int(v) for v in exact.as_integer_ratio()))
    if ex == 0:
        return Fraction(0) if y == 0 else Fraction(10**9)
    _, e = mpmath.frexp(exact)  # |exact| in [2^(e-1), 2^e)
    ulp = Fraction(2) ** (max(int(e) - 1, -1022) - 52)
    return abs(Fraction(y) - ex) / ulp


def exact_value(name: str, arg: tuple[float, ...], prec: int) -> mpf:
    with mp.workprec(prec):
        return MPMATH[name](*(mpf(a) for a in arg))


def main() -> None:
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 100_000
    cases = inputs(random.Random(SEED), n)
    results = evaluate(cases)
    print(f"Error of det against correctly rounded values, {n} inputs per function")
    print(f"(mpmath {mpmath.__version__}, seed {SEED}).")
    print()
    print("| Function | Max error | Correctly rounded | Worst input |")
    print("| --- | --- | --- | --- |")
    for name, args in cases.items():
        errors = []
        for arg, y in zip(args, results[name]):
            if not math.isfinite(y):
                raise SystemExit(f"{name}{arg} = {y}: non-finite result in the sample")
            errors.append((ulp_error(y, exact_value(name, arg, WORK_PREC_BITS)), arg, y))
        errors.sort(key=lambda t: t[0], reverse=True)
        for err, arg, y in errors[:CHECK_WORST]:
            check = ulp_error(y, exact_value(name, arg, CHECK_PREC_BITS))
            if abs(check - err) > Fraction(1, 10**9):
                raise SystemExit(f"{name}{arg}: precision disagreement {float(err)} vs {float(check)}")
        rounded = sum(1 for err, _, _ in errors if err <= Fraction(1, 2))
        worst_err, worst_arg, _ = errors[0]
        args_text = ", ".join(repr(a) for a in worst_arg)
        print(
            f"| {name} | {float(worst_err):.3f} ULP | {100 * rounded / len(errors):.2f} % "
            f"| `{name}({args_text})` |"
        )


if __name__ == "__main__":
    main()
