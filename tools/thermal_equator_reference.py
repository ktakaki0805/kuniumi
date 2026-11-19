#!/usr/bin/env python3
"""Reference values for the thermal equator, `kuniumi_sim::stages::climate::circulation` (ticket M4-21).

The thermal equator of a month is the warm centre of the zonal fit `f` of its
sea-level temperatures: the area-weighted mean latitude of the latitudes
within the cap (|phi| <= 60 deg) whose fit is within the band dT = 2 K of its
largest value f*, each weighted by its excess e = f - (f* - dT). With
z = sin(phi) (so that cos(phi) dphi = dz) that is

    phi_T = integral(asin(z) e(z) dz) / integral(e(z) dz)

over the intervals of z in [-sin 60 deg, sin 60 deg] where e > 0. This script
evaluates it for the zonally symmetric temperatures of the tests of
`crates/kuniumi-sim/src/stages/climate/circulation/tests.rs`, polynomials in z:
the largest value from the real roots of the derivative and the ends of the
cap, the ends of the band as real roots of f - (f* - dT), and the integrals by
mpmath's quadrature on each interval, at 50 significant digits. It is
independent of the code under test, which sums over latitudes 0.1 deg apart a
fit of the temperatures of a grid's cells.

Printed, in degrees:

- single peaks: parabolas in z, 300 - c (z - sin phi0)^2, and the fit warmest
  at the cap, 300 + 20 z;
- two competing peaks, 300 - 40 (z^2 - 1/4)^2 + delta z: peaks near +-30 deg,
  2.5 K above the valley at the equator, the northern warmer by about delta;
  the centre at a few delta, and the largest rise of the centre between two
  values of delta 0.02 K apart over -4 to 4 K.

Run from the repository root (see tools/README.md):

    tools/.venv/bin/python tools/thermal_equator_reference.py
"""

from __future__ import annotations

from mpmath import asin, degrees, mp, mpf, pi, polyroots, quad, sin

mp.dps = 50

CAP = pi / 3
BAND_K = mpf(2)
Z_CAP = sin(CAP)


def value(coefficients: list[mpf], z: mpf) -> mpf:
    """The polynomial with `coefficients` (highest degree first) at `z`."""
    total = mpf(0)
    for c in coefficients:
        total = total * z + c
    return total


def derivative(coefficients: list[mpf]) -> list[mpf]:
    n = len(coefficients) - 1
    return [c * (n - k) for k, c in enumerate(coefficients[:-1])]


def real_roots_in_cap(coefficients: list[mpf]) -> list[mpf]:
    # Leading zeros would trouble polyroots.
    while coefficients and coefficients[0] == 0:
        coefficients = coefficients[1:]
    if len(coefficients) < 2:
        return []
    roots = polyroots(coefficients, maxsteps=400, extraprec=400)
    return sorted(
        r.real for r in roots if abs(r.imag) < mpf(10) ** -30 and -Z_CAP < r.real < Z_CAP
    )


def centre_deg(coefficients: list[mpf]) -> tuple[mpf, mpf]:
    """The warm centre and the warmest latitude (degrees) of the profile."""
    coefficients = [mpf(c) for c in coefficients]
    candidates = [-Z_CAP, Z_CAP] + real_roots_in_cap(derivative(coefficients))
    z_top = max(candidates, key=lambda z: value(coefficients, z))
    top = value(coefficients, z_top)
    shifted = coefficients[:-1] + [coefficients[-1] - top + BAND_K]
    ends = [-Z_CAP] + real_roots_in_cap(shifted) + [Z_CAP]

    def excess(z: mpf) -> mpf:
        return value(shifted, z)

    numerator = mpf(0)
    denominator = mpf(0)
    for a, b in zip(ends, ends[1:]):
        if excess((a + b) / 2) > 0:
            numerator += quad(lambda z: asin(z) * excess(z), [a, b])
            denominator += quad(excess, [a, b])
    return degrees(numerator / denominator), degrees(asin(z_top))


def parabola(curvature: int, peak_deg: str) -> list[mpf]:
    """300 - c (z - s)^2 with s = sin(peak)."""
    s = sin(mpf(peak_deg) * pi / 180)
    c = mpf(curvature)
    return [-c, 2 * c * s, 300 - c * s * s]


def two_peaks(delta: mpf) -> list[mpf]:
    """300 - 40 (z^2 - 1/4)^2 + delta z = -40 z^4 + 20 z^2 + delta z + 297.5."""
    return [mpf(-40), mpf(0), mpf(20), mpf(delta), mpf("297.5")]


def main() -> None:
    print("single peaks: centre, warmest latitude (deg)")
    for name, profile in [
        ("300 - 40 (z - sin 12.3)^2", parabola(40, "12.3")),
        ("300 - 40 (z - sin -7)^2", parabola(40, "-7")),
        ("300 - 40 (z - sin 59)^2", parabola(40, "59")),
        ("300 - 30 (z - sin 20)^2", parabola(30, "20")),
        ("300 - 30 z^2", parabola(30, "0")),
        ("300 + 20 z", [mpf(20), mpf(300)]),
    ]:
        centre, warmest = centre_deg(profile)
        print(f"  {name}: {mp.nstr(centre, 12)} {mp.nstr(warmest, 12)}")
    print("two peaks, 300 - 40 (z^2 - 1/4)^2 + delta z: centre, warmest latitude (deg)")
    for delta in ["0", "0.5", "1", "2", "2.5", "4"]:
        for sign in ["", "-"] if delta != "0" else [""]:
            centre, warmest = centre_deg(two_peaks(mpf(sign + delta)))
            print(f"  delta {sign}{delta} K: {mp.nstr(centre, 12)} {mp.nstr(warmest, 12)}")
    largest, at = mpf(0), None
    previous = None
    for k in range(-200, 201):
        delta = mpf(k) / 50
        centre, _ = centre_deg(two_peaks(delta))
        if previous is not None and centre - previous > largest:
            largest, at = centre - previous, delta
        previous = centre
    print(
        f"  largest rise in a step of 0.02 K over -4..4 K: {mp.nstr(largest, 6)} deg, "
        f"up to delta {mp.nstr(at, 4)} K"
    )


if __name__ == "__main__":
    main()
