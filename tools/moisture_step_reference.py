#!/usr/bin/env python3
"""Reference values for the moisture's balance step, `kuniumi_sim::stages::climate::moisture` (ticket M4-22).

A still ocean planet: every cell at a uniform surface temperature T, no wind
(the bulk formula's wind speed is the gust speed U_g alone, no transport, no
lift), and the precipitation threshold r_c = 1, so that a cell below its
saturated column W_s never rains. Each cell is then alone and linear,

    dW/dt = k (W_s - W),   k = C_E U_g / H_w,

with W_s = rho_vs H_w, and its periodic state is W_s. The model steps it by
backward Euler, 48 sub-steps a year of length dt = P/48,

    W' = (W/dt + k W_s) / (1/dt + k),

and above the threshold (only after a step overshoots W_s) with the excess
raining in tau_c, W' = (W/dt + k W_s + W_s/tau_c) / (1/dt + k + 1/tau_c).
A year is periodic when no monthly mean (of four sub-steps) changes by more
than 1e-6 kg/m2 from the year before, at most 10 years. From W = 0 the years
alone need about (ln(W_s / 1e-6)) / (kP) years; with the balance step after
each year that is not periodic, from x at the start of the year to y at its
end, a cell alone has

    Y/a = dt sum_s (k + [raining in s] / tau_c) + mu      (mu = 1e-9)
    u = (y - x) / (Y/a)                                  (the mean shift)
    v_0 = u,  v_{j+1} = v_j / (1 + (Y/a)/48),  j < 24     (half a year)
    next start max(y + v_24, 0).

This script runs that recurrence at 50 significant digits with the closed
forms of the saturation (Clausius-Clapeyron from the triple point with the
latent heat of vaporisation, the scale height of a saturated column along
the lapse rate Gamma), independent of the code under test, which solves the
same rows of every cell of a grid by its sweeps, its active set method and
its sparse solves. Printed: W_s, k, kP, the years and the last change with
and without the step, and the twelve monthly means of the last year with it.

It also prints the factor g(z) by which a year and its step leave an error
of the start of a linear cell alone, z = kP (without the mass mu): the years
leave mu_y = (1 + z/48)^-48, and the step v_24 = u (1 + z/48)^-24 with
u = (y - x)/z, so g = mu_y - (1 - mu_y)(1 + z/48)^-24 / z; with the weight
1/(1 + z/2) instead of the half year (Pade), g = mu_y - (1 - mu_y)/(z (1 +
z/2)); with the mean shift u alone, g = mu_y - (1 - mu_y)/z. Printed: the
largest |g| over z from 10^-3 to 10^4 (1,401 values, logarithmically spaced)
and g at z = 0.155 (corner 93's slow cell, ticket M4-22).

Run from the repository root (see tools/README.md):

    tools/.venv/bin/python tools/moisture_step_reference.py
"""

from __future__ import annotations

from mpmath import exp, mp, mpf, nstr

mp.dps = 50

# The model's constants (moisture.rs, energy_balance.rs).
TRIPLE_POINT_K = mpf("273.16")
TRIPLE_POINT_PRESSURE_PA = mpf("611.657")
LATENT_HEAT_J_KG = mpf("2.501e6")
VAPOUR_GAS_CONSTANT = mpf("461.5")
LAPSE_RATE_K_M = mpf("6.5e-3")
DALTON_NUMBER = mpf("1.2e-3")
GUST_SPEED_M_S = mpf(2)
CONDENSATION_TIME_S = mpf(14400)
TOLERANCE_KG_M2 = mpf("1e-6")
MAX_YEARS = 10
BALANCE_MASS = mpf("1e-9")
MONTHS = 12
SUBSTEPS_PER_MONTH = 4
SUBSTEPS = MONTHS * SUBSTEPS_PER_MONTH

# The planet of the test.
SURFACE_K = mpf(288)
YEAR_S = mpf(200000)


def saturated_column() -> tuple[mpf, mpf]:
    """W_s (kg/m2) and H_w (m) at SURFACE_K."""
    t = SURFACE_K
    l, r = LATENT_HEAT_J_KG, VAPOUR_GAS_CONSTANT
    e = TRIPLE_POINT_PRESSURE_PA * exp(l / r * (1 / TRIPLE_POINT_K - 1 / t))
    rho = e / (r * t)
    h = r * t * t / (LAPSE_RATE_K_M * (l - r * t))
    return rho * h, h


def run(step: bool, max_years: int) -> tuple[int, mpf, list[mpf]]:
    """The periodic state from W = 0: years, last change, monthly means."""
    ws, h = saturated_column()
    k = DALTON_NUMBER * GUST_SPEED_M_S / h
    dt = YEAR_S / SUBSTEPS
    w = mpf(0)
    previous: list[mpf] | None = None
    years = 0
    while True:
        years += 1
        x = w
        wet = 0
        means = []
        for _ in range(MONTHS):
            total = mpf(0)
            for _ in range(SUBSTEPS_PER_MONTH):
                dry = (w / dt + k * ws) / (1 / dt + k)
                if dry > ws:
                    w = (w / dt + k * ws + ws / CONDENSATION_TIME_S) / (
                        1 / dt + k + 1 / CONDENSATION_TIME_S
                    )
                    wet += 1
                else:
                    w = dry
                total += w
            means.append(total / SUBSTEPS_PER_MONTH)
        if previous is None:
            change = mp.inf
        else:
            change = max(abs(a - b) for a, b in zip(means, previous))
        if change <= TOLERANCE_KG_M2 or years >= max_years:
            return years, change, means
        previous = means
        if step:
            z = dt * (SUBSTEPS * k + wet / CONDENSATION_TIME_S) + BALANCE_MASS
            v = (w - x) / z
            for _ in range(SUBSTEPS // 2):
                v = v / (1 + z / SUBSTEPS)
            w = max(w + v, mpf(0))


def contraction(z: mpf, kind: str) -> mpf:
    """The factor g(z) of an error of a cell alone after a year and its step."""
    years = (1 + z / SUBSTEPS) ** (-SUBSTEPS)
    if kind == "the plan":
        weight = (1 + z / SUBSTEPS) ** (-(SUBSTEPS // 2)) / z
    elif kind == "Pade":
        weight = 1 / (z * (1 + z / 2))
    else:
        weight = 1 / z
    return years - (1 - years) * weight


def main() -> None:
    ws, h = saturated_column()
    k = DALTON_NUMBER * GUST_SPEED_M_S / h
    print(f"W_s = {nstr(ws, 17)} kg/m2, H_w = {nstr(h, 17)} m")
    print(f"k = {nstr(k, 17)} 1/s, kP = {nstr(k * YEAR_S, 17)}")
    years, change, _ = run(False, 1000)
    print(f"without the step: {years} years, last change {nstr(change, 6)} kg/m2")
    years, change, _ = run(False, MAX_YEARS)
    print(f"without the step, {MAX_YEARS} years: last change {nstr(change, 6)} kg/m2")
    years, change, means = run(True, MAX_YEARS)
    print(f"with the step: {years} years, last change {nstr(change, 17)} kg/m2")
    for m, value in enumerate(means):
        print(f"  month {m + 1}: {nstr(value, 20)} kg/m2 (W_s - it: {nstr(ws - value, 6)})")
    zs = [mpf(10) ** (mpf(j) / 200 - 3) for j in range(1401)]
    for kind in ("the plan", "Pade", "the mean shift"):
        worst = max(zs, key=lambda z: abs(contraction(z, kind)))
        print(
            f"{kind}: largest |g| {nstr(abs(contraction(worst, kind)), 4)} at z = "
            f"{nstr(worst, 4)}; g(0.155) = {nstr(contraction(mpf('0.155'), kind), 4)}"
        )


if __name__ == "__main__":
    main()
