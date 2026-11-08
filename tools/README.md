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

Regenerate and compare with the committed file:

```sh
tools/.venv/bin/python tools/det_reference.py | diff - crates/kuniumi-math/src/det/golden.rs
```

Measure the error of `det` against the exact value (10⁶ inputs per function
take about three minutes; `det_error.py` runs cargo, so call it through mise):

```sh
mise exec -- tools/.venv/bin/python tools/det_error.py 1000000
```
