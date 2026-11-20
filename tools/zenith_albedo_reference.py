#!/usr/bin/env python3
"""Reference values of the albedo by the solar zenith angle (ticket M4-10).

Prints the numbers that the tests of `kuniumi_sim::stages::climate::insolation`
(the insolation-weighted cosine of the zenith angle) and of
`kuniumi_sim::stages::climate::energy_balance` (the zenith-angle law of the
ice-free albedo and its calibration) compare with, and the constants of the
law, computed with mpmath (30 digits):

- the insolation-weighted mean cosine of the zenith angle over a day,
  `integral of cos^2 z dh / integral of cos z dh` over the hour angles of
  daylight, `cos z = sin phi sin delta + cos phi cos delta cos h`, integrated
  numerically (Gauss-Legendre), not with the closed forms of the code;
- the monthly value, the same ratio of time integrals over the twelve equal
  months of the calendar, weighted by the flux `(a / r)^2`, with the orbit of
  `insolation_reference.py` (Kepler's equation solved at each instant, the
  months split where polar day or night begins or ends);
- the constants `b` and `mu_r` of the law `alpha_free = alpha_type + b (mu_r - mu)`:
  the annual ice-free absorbed sunlight of Earth's obliquity has the global
  mean and the P2 moment of the form of ticket M4-02,
  `alpha_type + 0.08 P2(sin phi)`. The annual integrals over the year do not
  depend on the eccentricity, the perihelion or the calendar (Kepler's second
  law: `(a/r)^2 dt` is proportional to the star's longitude), so they are
  integrals over the longitude of a circular orbit; the day integrals here are
  the closed forms `F = h0 A + B sin h0` and
  `G = A^2 h0 + 2 A B sin h0 + B^2 (h0 + sin h0 cos h0) / 2`, checked against
  the quadrature over the hour angle first;
- with the rounded constants of the code (0.338, 0.644): the two moments
  against those of the P2 form, the insolation-weighted annual ice-free
  albedo of each 10 degree band (relative to `alpha_type`) for both, and at
  obliquity 90 degrees the annual weighted cosine and the ocean's ice-free
  co-albedo at the pole and the equator (hand values pi/4 and pi^2/16).

Run from the repository root (see tools/README.md); it takes one to three
minutes:

    tools/.venv/bin/python tools/zenith_albedo_reference.py
"""

from __future__ import annotations

import os
import sys

from mpmath import mp, mpf

# The orbit of the insolation's reference script, imported from beside this
# one without leaving a __pycache__ in tools/.
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from insolation_reference import Orbit  # noqa: E402

mp.dps = 30

ALPHA_P2 = mpf("0.08")
OCEAN_ALBEDO = mpf("0.263")
SLOPE = mpf("0.338")
REFERENCE = mpf("0.644")


def daylight(a: mpf, b: mpf) -> mpf:
    """The hour angle of sunset for cos z = a + b cos h (pi in polar day, 0 in polar night)."""
    if a >= b:
        return mp.pi
    if -a >= b:
        return mpf(0)
    return mp.acos(-a / b)


def day_moments(phi: mpf, delta: mpf) -> tuple[mpf, mpf]:
    """Integrals of cos z and cos^2 z over the daylight hour angles, by quadrature."""
    a = mp.sin(phi) * mp.sin(delta)
    b = mp.cos(phi) * mp.cos(delta)
    h0 = daylight(a, b)
    if h0 == 0:
        return mpf(0), mpf(0)
    f = mp.quad(lambda h: a + b * mp.cos(h), [-h0, h0], method="gauss-legendre")
    g = mp.quad(lambda h: (a + b * mp.cos(h)) ** 2, [-h0, h0], method="gauss-legendre")
    return f, g


def daily_cosine(phi_deg: str, delta_deg: str) -> mpf:
    f, g = day_moments(mp.radians(mpf(phi_deg)), mp.radians(mpf(delta_deg)))
    return g / f if f > 0 else mpf(0)


def closed_forms(a: mpf, b: mpf) -> tuple[mpf, mpf]:
    """The half-day integrals F and G in closed form."""
    h0 = daylight(a, b)
    if h0 == 0:
        return mpf(0), mpf(0)
    s, c = mp.sin(h0), mp.cos(h0)
    return h0 * a + b * s, a * a * h0 + 2 * a * b * s + b * b * (h0 + s * c) / 2


def monthly_cosine(orbit: Orbit, phi: mpf, month: int) -> mpf:
    """The insolation-weighted mean cosine of month `month` (1 to 12)."""
    t0, t1 = mpf(month - 1) / 12, mpf(month) / 12
    points = [t0] + orbit.kinks(phi, t0, t1) + [t1]

    def moments(t: mpf) -> tuple[mpf, mpf]:
        sd, rho2 = orbit.at(t)
        f, g = day_moments(phi, mp.asin(sd))
        return rho2 * f, rho2 * g

    total_f = total_g = mpf(0)
    for a, b in zip(points, points[1:]):
        total_f += mp.quad(lambda t: moments(t)[0], [a, b], method="gauss-legendre")
        total_g += mp.quad(lambda t: moments(t)[1], [a, b], method="gauss-legendre")
    return total_g / total_f if total_f > 0 else mpf(0)


class Annual:
    """Annual integrals over the star's longitude of F and G at sin(latitude) x."""

    def __init__(self, obliquity_deg: str):
        self.eps = mp.radians(mpf(obliquity_deg))
        self.cache: dict = {}

    def __call__(self, x: mpf) -> tuple[mpf, mpf]:
        if x not in self.cache:
            cl = mp.sqrt(1 - x * x)
            se = mp.sin(self.eps)

            def at(lam: mpf, which: int) -> mpf:
                sd = se * mp.sin(lam)
                return closed_forms(x * sd, cl * mp.sqrt(1 - sd * sd))[which]

            # Split at the solstices, where the declination turns (and at the
            # equator of obliquity 90 degrees cos(delta) = |cos(lambda)| has a
            # kink), and where polar day or night begins or ends.
            points = [k * mp.pi / 2 for k in range(5)]
            if cl < se:
                l0 = mp.asin(cl / se)
                points += [l0, mp.pi - l0, mp.pi + l0, 2 * mp.pi - l0]
            points = sorted(points)
            self.cache[x] = (
                mp.quad(lambda lam: at(lam, 0), points),
                mp.quad(lambda lam: at(lam, 1), points),
            )
        return self.cache[x]

    def integral(self, fun, lo: mpf, hi: mpf) -> mpf:
        """Integral over x in [lo, hi] (0 <= lo < hi <= 1), split at the polar circle."""
        pc = mp.cos(self.eps)
        points = [lo] + ([pc] if lo < pc < hi else []) + [hi]
        return mp.quad(fun, points)


def p2(x: mpf) -> mpf:
    return (3 * x * x - 1) / 2


def show(name: str, value: mpf) -> None:
    print(f"{name} = {mp.nstr(value, 17)}")


def main() -> None:
    print("# Closed forms of the day integrals against the quadrature (largest difference).")
    worst = mpf(0)
    for phi, delta in [("0", "0"), ("45", "10"), ("70", "20"), ("-30", "20"), ("89", "1.5"), ("60", "-29.5"), ("90", "23.44")]:
        p, d = mp.radians(mpf(phi)), mp.radians(mpf(delta))
        a, b = mp.sin(p) * mp.sin(d), mp.cos(p) * mp.cos(d)
        f, g = closed_forms(a, b)
        fq, gq = day_moments(p, d)
        worst = max(worst, abs(2 * f - fq), abs(2 * g - gq))
    show("closed_forms_error", worst)

    print("# Daily insolation-weighted cosines: latitude, declination (degrees).")
    for phi, delta in [("0", "0"), ("45", "10"), ("70", "20"), ("-30", "20"), ("89", "1.5"), ("60", "-29.5"), ("-80", "-15"), ("90", "23.44"), ("23.44", "23.44")]:
        show(f"daily({phi}, {delta})", daily_cosine(phi, delta))

    earth = Orbit("1361", "23.44", "0.0167", "102.94")
    eccentric = Orbit("1361", "60", "0.5", "250")
    for name, orbit in [("earth", earth), ("eccentric", eccentric)]:
        print(f"# {name}: monthly insolation-weighted cosines, latitude (degrees), month.")
        for phi in ["90", "60", "0", "-45", "-80"]:
            for month in [1, 4, 6, 7, 10, 12]:
                show(f"{name}.monthly({phi}, {month})", monthly_cosine(orbit, mp.radians(mpf(phi)), month))

    print("# The law's constants from Earth's obliquity (23.44 degrees).")
    annual = Annual("23.44")
    zero, one = mpf(0), mpf(1)
    f0 = annual.integral(lambda x: annual(x)[0], zero, one)
    f2 = annual.integral(lambda x: p2(x) * annual(x)[0], zero, one)
    h = annual.integral(lambda x: p2(x) ** 2 * annual(x)[0], zero, one)
    g0 = annual.integral(lambda x: annual(x)[1], zero, one)
    g2 = annual.integral(lambda x: p2(x) * annual(x)[1], zero, one)
    # u = b mu_r: u F0 - b G0 = a2 F2 and u F2 - b G2 = a2 H.
    det = -f0 * g2 + g0 * f2
    u = ALPHA_P2 * (-f2 * g2 + g0 * h) / det
    slope = ALPHA_P2 * (f0 * h - f2 * f2) / det
    show("b", slope)
    show("mu_r", u / slope)
    show("Z_r_deg", mp.degrees(mp.acos(u / slope)))

    print("# With b = 0.338 and mu_r = 0.644: the moments of the annual ice-free albedo")
    print("# times the insolation, against the P2 form's (relative difference).")
    new0 = SLOPE * (REFERENCE * f0 - g0)
    new2 = SLOPE * (REFERENCE * f2 - g2)
    old0 = ALPHA_P2 * f2
    old2 = ALPHA_P2 * h
    # The P2 form's global moment of absorbed sunlight is a2 F2 relative to
    # the whole absorbed sunlight (1 - alpha_type) F0; compare the absorbed
    # sunlight of an ice-free ocean planet.
    show("global_moment_relative", ((1 - OCEAN_ALBEDO) * f0 - new0) / ((1 - OCEAN_ALBEDO) * f0 - old0) - 1)
    show("p2_moment_relative", ((1 - OCEAN_ALBEDO) * f2 - new2) / ((1 - OCEAN_ALBEDO) * f2 - old2) - 1)

    print("# Insolation-weighted annual ice-free albedo less alpha_type, 10 degree bands (either hemisphere).")
    for lo in range(0, 90, 10):
        x0, x1 = mp.sin(mp.radians(lo)), mp.sin(mp.radians(lo + 10))
        q = annual.integral(lambda x: annual(x)[0], x0, x1)
        new = annual.integral(lambda x: SLOPE * (REFERENCE * annual(x)[0] - annual(x)[1]), x0, x1) / q
        old = annual.integral(lambda x: ALPHA_P2 * p2(x) * annual(x)[0], x0, x1) / q
        show(f"band({lo}).zenith", new)
        show(f"band({lo}).p2", old)

    print("# Obliquity 90 degrees: annual weighted cosine and the ocean's ice-free co-albedo.")
    oblique = Annual("90")
    for name, x, hand in [("pole", mpf(1), mp.pi / 4), ("equator", mpf(0), mp.pi**2 / 16)]:
        f, g = oblique(x)
        mu = g / f
        show(f"90.{name}.mu", mu)
        show(f"90.{name}.mu_hand_difference", mu - hand)
        show(f"90.{name}.co_albedo", 1 - OCEAN_ALBEDO - SLOPE * (REFERENCE - mu))
        show(f"90.{name}.co_albedo_p2_form", 1 - OCEAN_ALBEDO - ALPHA_P2 * p2(x))


if __name__ == "__main__":
    main()
