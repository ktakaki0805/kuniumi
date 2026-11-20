# Development tools

Python scripts for development only (D-11): computing independent reference
values for tests and, later, converting reference data. They are never part of
the build, and the Rust code never runs them.

Python is pinned in `mise.toml` and the packages in `requirements.txt`. Set up
a virtual environment once, from the repository root:

```sh
mise install
mise exec -- python -m venv tools/.venv
tools/.venv/bin/pip install -r tools/requirements.txt
```

## Scripts

| Script | Output | Ticket |
| --- | --- | --- |
| `det_reference.py` | `crates/kuniumi-math/src/det/golden.rs`: correctly rounded values of the `det` functions (mpmath) | M1-01 |
| `det_error.py` | Table of the error of the `det` functions against the exact value (mpmath), quoted in the `det` documentation | M1-01 |
| `crust_reference.py` | The age at which the two depth curves of Parsons & Sclater (1977) meet, depths at a few ages and heights by Airy isostasy (mpmath), quoted in the tests of `kuniumi_sim::stages::crust` | M2-05 |
| `hypsometry_reference.py` | The numbers of `metrics.txt` of a world directory (land fraction, means, extremes, percentiles, modes and valley of the elevation histogram) with exact rational arithmetic, quoted in the tests of `crates/kuniumi-cli/tests/cli.rs`; standard library only | M2-06 |
| `colormap_tables.py` | The viridis, RdBu and Paired tables of `kuniumi_io::colormap` from their sources (downloaded at pinned commits, SHA-256 checked), the viridis reference colours of its tests and the smallest colour difference in Paired; standard library only | M1-08 |
| `insolation_reference.py` | Daily, monthly and annual means of the insolation for Earth and an eccentric orbit, and the dates of Earth's perihelion, equinoxes and solstices, by integration in time with Kepler's equation and over the hour angle (mpmath), quoted in the tests of `kuniumi_sim::stages::climate::insolation` | M4-01 |
| `thermal_equator_reference.py` | The warm centre of the thermal equator (the centre of the latitudes within 2 K of the warmest, ±60°) of zonally symmetric temperatures: single peaks and two competing peaks, with the largest step of the centre (mpmath), quoted in the tests of `kuniumi_sim::stages::climate::circulation` | M4-21 |
| `moisture_step_reference.py` | The periodic state of a still ocean planet whose cells never rain (each alone and linear) by the recurrence of the moisture's sub-steps and its balance step, with the years, the last change and the monthly vapour, and the factor by which a year and its step leave an error of a cell alone (mpmath), quoted in the tests and the documentation of `kuniumi_sim::stages::climate::moisture` | M4-22 |
| `zenith_albedo_reference.py` | The insolation-weighted cosine of the zenith angle by day (quadrature over the hour angle) and by month (time integrals with Kepler's equation, the orbit of `insolation_reference.py`), the constants `b` and `μ_r` of the zenith angle law of the albedo from the moments of North & Coakley's `P₂` form on Earth's obliquity, the bands of both and the values at obliquity 90° (mpmath), quoted in the tests and the documentation of `kuniumi_sim::stages::climate::insolation` and `energy_balance` | M4-10 |
| `font8x8_table.py` | The glyph tables of the viewer's bitmap font (`crates/kuniumi-viewer/src/font.rs`) from font8x8 (downloaded at a pinned commit, SHA-256 checked), and the hash of the whole table that its test quotes; standard library only | MV-06 |

Regenerate and compare with the committed file:

```sh
tools/.venv/bin/python tools/det_reference.py | diff - crates/kuniumi-math/src/det/golden.rs
```

Measure the error of `det` against the exact value (10⁶ inputs per function
take about three minutes; `det_error.py` runs cargo, so call it through mise):

```sh
mise exec -- tools/.venv/bin/python tools/det_error.py 1000000
```

Print the colour tables of `kuniumi_io::colormap` from their sources (needs
network access; the standard library is enough) and compare them with
`crates/kuniumi-io/src/colormap.rs` and its tests:

```sh
mise exec -- python tools/colormap_tables.py
```

Print the tables of the viewer's bitmap font from their source (needs
network access) and compare them with `font.rs`:

```sh
mise exec -- python tools/font8x8_table.py | head -n 98 | diff - <(sed -n '/^const BASIC/,/^const DEGREE/p' crates/kuniumi-viewer/src/font.rs | grep -v '^$\|^///')
```

Print the reference values of the Crust stage's tests:

```sh
tools/.venv/bin/python tools/crust_reference.py
```

Print the reference values of the insolation's tests (about 15 s):

```sh
tools/.venv/bin/python tools/insolation_reference.py
```

Print the reference values of the thermal equator's tests (about 5 s):

```sh
tools/.venv/bin/python tools/thermal_equator_reference.py
```

Print the reference values and the constants of the zenith angle law of the
albedo (one to three minutes):

```sh
tools/.venv/bin/python tools/zenith_albedo_reference.py
```

Print the reference metrics of the golden world and compare them with its
`metrics.txt` (the standard library is enough):

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 10000 --out out/M2-06/10k
mise exec -- python tools/hypsometry_reference.py out/M2-06/10k
```
