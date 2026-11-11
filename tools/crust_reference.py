#!/usr/bin/env python3
"""Reference values for the Crust stage, `kuniumi_sim::stages::crust` (ticket M2-05).

Prints, at 100 significant digits (mpmath, 400 bits), the values that the
tests of `crates/kuniumi-sim/src/stages/crust/tests.rs` compare with:

- the age t* (Myr) at which the two depth curves of Parsons & Sclater (1977)
  meet, `2500 + 350 sqrt(t) = 6400 - 3200 exp(-t / 62.8)` (m), found by
  bisection in [20, 70] Myr, where the first is below the second at 20 Myr and
  above it at 70 Myr; and the depth d* = d(t*) there;
- the depth of the curve switched at t* at a few ages;
- the height by Airy isostasy, `(T - 35 km) (3300 - 2800) / 3300 + 300 m`,
  of crust 15, 35 and 70 km thick, as exact fractions too.

Run from the repository root (see tools/README.md):

    tools/.venv/bin/python tools/crust_reference.py
"""

from __future__ import annotations

from fractions import Fraction

from mpmath import mp, mpf

mp.prec = 400


def young(t: mpf) -> mpf:
    return 2500 + 350 * mp.sqrt(t)


def old(t: mpf) -> mpf:
    return 6400 - 3200 * mp.exp(-t / mpf("62.8"))


def switch() -> mpf:
    lo, hi = mpf(20), mpf(70)
    assert young(lo) < old(lo) and young(hi) > old(hi)
    for _ in range(600):
        mid = (lo + hi) / 2
        if young(mid) < old(mid):
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def depth(t: mpf, t_switch: mpf) -> mpf:
    return young(t) if t < t_switch else old(t)


def airy(thickness_m: int) -> Fraction:
    return Fraction(thickness_m - 35_000) * Fraction(3300 - 2800, 3300) + 300


def main() -> None:
    t_switch = switch()
    print(f"switch_myr = {mp.nstr(t_switch, 100)}")
    print(f"switch_depth_m = {mp.nstr(young(t_switch), 100)}")
    print(f"young - old at the switch = {mp.nstr(young(t_switch) - old(t_switch), 5)}")
    for t in ["0", "1", "4", "16", "25", "30", "60", "100", "180"]:
        print(f"depth_m({t} Myr) = {mp.nstr(depth(mpf(t), t_switch), 100)}")
    for thickness_m in [15_000, 35_000, 70_000]:
        h = airy(thickness_m)
        print(f"airy_m({thickness_m} m) = {h} = {mp.nstr(mpf(h.numerator) / h.denominator, 100)}")


if __name__ == "__main__":
    main()
