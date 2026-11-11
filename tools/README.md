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
| `colormap_tables.py` | The viridis, RdBu and Paired tables of `kuniumi_io::colormap` from their sources (downloaded at pinned commits, SHA-256 checked), the viridis reference colours of its tests and the smallest colour difference in Paired; standard library only | M1-08 |

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

Print the reference values of the Crust stage's tests:

```sh
tools/.venv/bin/python tools/crust_reference.py
```
