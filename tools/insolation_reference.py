#!/usr/bin/env python3
"""Reference values of the insolation, `kuniumi_sim::stages::climate::insolation` (ticket M4-01).

Prints the numbers that `crates/kuniumi-sim/src/stages/climate/insolation/tests.rs`
compares with, computed with mpmath (30 digits) along a route independent of
the Rust code, which integrates the daily mean over the star's longitude with a
closed form:

- the daily mean at latitude phi and declination delta is the mean over one
  turn of the hour angle h of the flux on a horizontal surface,
  `S rho2 / (2 pi) * integral of max(0, sin phi sin delta + cos phi cos delta cos h) dh`,
  integrated numerically (Gauss-Legendre) between the hour angles of sunrise
  and sunset, not with Berger's closed form; `rho2` is `(a / r)^2`;
- the position on the orbit at a year fraction tau comes from the time since
  periapsis: the year fraction of the northern spring equinox is fixed at
  79.5 / 365 (PMIP), the time from periapsis to the equinox is the integral
  of `dt / dv = (1 - e^2)^(3/2) / (2 pi (1 + e cos v)^2)` (years per radian
  of true anomaly) up to the equinox's true anomaly `v = pi - varpi`, and at
  any tau Kepler's equation is solved with `mp.findroot` (Illinois method on
  the bracket [M - e, M + e]); then `r / a = 1 - e cos E`,
  `v = atan2(sqrt(1 - e^2) sin E, cos E - e)` and the star's longitude
  `lambda = v + varpi + pi`, `sin delta = sin(eps) sin(lambda)`;
- monthly means are time integrals over the twelve equal months (year
  fractions k / 12 to (k + 1) / 12), split where polar day or night begins
  or ends (sign changes on 64 points per month, then bisection on tau), so
  that every piece is smooth; annual means are the means of the twelve.

Also: the year fractions of periapsis, of the equinoxes and of the solstices
(from the time-since-periapsis integral at the true anomaly of each), and the
annual means of the closed forms that the tests quote as hand values (the
pole `S sin(eps) / (pi sqrt(1 - e^2))`, the equator
`S / pi * 2 / pi * E(sin eps) / sqrt(1 - e^2)` with the complete elliptic
integral of the second kind).

Run from the repository root (see tools/README.md); it takes about 15 seconds:

    tools/.venv/bin/python tools/insolation_reference.py
"""

from __future__ import annotations

from mpmath import mp, mpf

mp.dps = 30

SPRING_EQUINOX = mpf("79.5") / 365


class Orbit:
    def __init__(self, flux: str, obliquity_deg: str, eccentricity: str, perihelion_deg: str):
        self.s = mpf(flux)
        self.eps = mp.radians(mpf(obliquity_deg))
        self.e = mpf(eccentricity)
        self.varpi = mp.radians(mpf(perihelion_deg))
        # Year fraction of periapsis: the equinox minus the time from
        # periapsis to the equinox.
        self.periapsis = SPRING_EQUINOX - self.time_from_periapsis(mp.pi - self.varpi)

    def time_from_periapsis(self, v: mpf) -> mpf:
        """Years from periapsis to true anomaly v, from Kepler's second law."""
        e = self.e
        return mp.quad(lambda u: (1 - e * e) ** mpf(1.5) / (2 * mp.pi * (1 + e * mp.cos(u)) ** 2), [0, v])

    def eccentric_anomaly(self, m: mpf) -> mpf:
        e = self.e
        if e == 0:
            return m
        return mp.findroot(lambda x: x - e * mp.sin(x) - m, (m - e, m + e), solver="illinois")

    def at(self, tau: mpf) -> tuple[mpf, mpf]:
        """(sin delta, (a / r)^2) at year fraction tau."""
        e = self.e
        m = 2 * mp.pi * (tau - self.periapsis)
        big_e = self.eccentric_anomaly(m)
        r = 1 - e * mp.cos(big_e)
        v = mp.atan2(mp.sqrt(1 - e * e) * mp.sin(big_e), mp.cos(big_e) - e)
        lam = v + self.varpi + mp.pi
        return mp.sin(self.eps) * mp.sin(lam), 1 / (r * r)

    def daily(self, phi: mpf, tau: mpf) -> mpf:
        sd, rho2 = self.at(tau)
        return daily_mean(self.s * rho2, phi, mp.asin(sd))

    def mean(self, phi: mpf, t0: mpf, t1: mpf) -> mpf:
        """Time mean of the daily mean at phi over year fractions [t0, t1]."""
        points = [t0] + self.kinks(phi, t0, t1) + [t1]
        total = mpf(0)
        for a, b in zip(points, points[1:]):
            total += mp.quad(lambda t: self.daily(phi, t), [a, b], method="gauss-legendre")
        return total / (t1 - t0)

    def kinks(self, phi: mpf, t0: mpf, t1: mpf) -> list[mpf]:
        """Year fractions in (t0, t1) where polar day (sin phi sin delta = cos phi cos delta)
        or polar night (sin phi sin delta = -cos phi cos delta) begins or ends at phi."""

        def boundary(sign: int):
            def g(t: mpf) -> mpf:
                sd, _ = self.at(t)
                return mp.sin(phi) * sd - sign * mp.cos(phi) * mp.sqrt(1 - sd * sd)

            return g

        out = []
        n = 64
        grid = [t0 + (t1 - t0) * k / n for k in range(n + 1)]
        for g in [boundary(1), boundary(-1)]:
            values = [g(t) for t in grid]
            for a, b, ga, gb in zip(grid, grid[1:], values, values[1:]):
                if ga * gb < 0:
                    lo, hi = a, b
                    for _ in range(120):
                        mid = (lo + hi) / 2
                        if g(lo) * g(mid) <= 0:
                            hi = mid
                        else:
                            lo = mid
                    out.append((lo + hi) / 2)
        return sorted(out)

    def fraction_at_longitude(self, lam: mpf) -> mpf:
        """Year fraction in [0, 1) at which the star's longitude is lam."""
        v = lam - self.varpi - mp.pi
        v = v - 2 * mp.pi * mp.floor(v / (2 * mp.pi))
        t = self.periapsis + self.time_from_periapsis(v)
        return t - mp.floor(t)


def daily_mean(flux: mpf, phi: mpf, delta: mpf) -> mpf:
    """Mean over a day of the flux on a horizontal surface, by integration over the hour angle."""
    a = mp.sin(phi) * mp.sin(delta)
    b = mp.cos(phi) * mp.cos(delta)
    if a >= b:
        h0 = mp.pi
    elif -a >= b:
        return mpf(0)
    else:
        h0 = mp.acos(-a / b)
    integral = mp.quad(lambda h: a + b * mp.cos(h), [-h0, h0], method="gauss-legendre")
    return flux * integral / (2 * mp.pi)


def show(name: str, value: mpf) -> None:
    print(f"{name} = {mp.nstr(value, 17)}")


def main() -> None:
    deg = lambda x: mp.radians(mpf(x))  # noqa: E731

    print("# Daily means (W/m2): flux S (a/r)^2, latitude, declination.")
    for flux, phi, delta in [
        ("1361", "0", "0"),
        ("1361", "90", "23.44"),
        ("1361", "-90", "-23.44"),
        ("1361", "90", "0"),
        ("1361", "45", "10"),
        ("1361", "-30", "20"),
        ("1361", "70", "-15"),
        ("1361", "70", "20"),
        ("1300", "-80", "-15"),
        ("2000", "60", "-29.5"),
        ("1361", "60", "-30"),
        ("400", "89", "1.5"),
    ]:
        show(f"daily({flux}, {phi}, {delta})", daily_mean(mpf(flux), deg(phi), deg(delta)))

    earth = Orbit("1361", "23.44", "0.0167", "102.94")
    eccentric = Orbit("1361", "60", "0.5", "250")
    for name, orbit in [("earth", earth), ("eccentric", eccentric)]:
        print(f"# {name}: year fractions of periapsis and of the star's longitude 0, 90, 180, 270 deg.")
        show(f"{name}.periapsis", orbit.periapsis - mp.floor(orbit.periapsis))
        for lam in ["0", "90", "180", "270"]:
            show(f"{name}.fraction(lambda = {lam})", orbit.fraction_at_longitude(deg(lam)))
        print(f"# {name}: annual means (W/m2).")
        for phi in ["0", "90", "-90", "60", "-45"]:
            months = [orbit.mean(deg(phi), mpf(k) / 12, mpf(k + 1) / 12) for k in range(12)]
            show(f"{name}.annual({phi})", sum(months) / 12)
        print(f"# {name}: monthly means (W/m2), months 1 to 12.")
        for phi in ["90", "-90", "60", "0", "-45"]:
            for month in [1, 3, 6, 7, 9, 12]:
                t0 = mpf(month - 1) / 12
                t1 = mpf(month) / 12
                show(f"{name}.monthly({phi}, {month})", orbit.mean(deg(phi), t0, t1))

    print("# Closed forms of the annual means at e = 0.0167 (W/m2).")
    e = mpf("0.0167")
    for eps in ["23.44", "45", "60", "90"]:
        k = mp.sin(deg(eps))
        show(f"pole({eps})", 1361 * k / (mp.pi * mp.sqrt(1 - e * e)))
        show(f"equator({eps})", 1361 / mp.pi * 2 / mp.pi * mp.ellipe(k * k) / mp.sqrt(1 - e * e))


if __name__ == "__main__":
    main()
