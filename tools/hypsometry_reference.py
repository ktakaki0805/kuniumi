#!/usr/bin/env python3
"""Reference values of the elevation metrics, `kuniumi_sim::metrics` (ticket M2-06).

Reads a world directory written by `kuniumi generate` (`surface_height.kfield`,
`cell_area.kfield` and the sea level in `world.json`) and prints, as `key =
value` lines, the numbers that `metrics.txt` of that directory reports, computed
independently of the Rust code:

- every elevation `h - z` and every bin `floor((h - z) / 100 m)` exactly, with
  rational arithmetic (`fractions.Fraction` of the doubles);
- the land fraction, the means and the shares of the bins from exact sums of the
  areas (`Fraction`), rounded once at the end;
- the percentiles from the cells sorted by their exact elevation (ties by cell
  index) and exact running sums of the areas, against `p` times the total;
- the modes (the bins with the largest share below and above the sea level, the
  lowest of equal ones), the valley between them (the smallest bin, the lowest of
  equal ones) and whether it holds at most half of the smaller mode.

A cell is ocean where its height lies below the sea level and land elsewhere.
The tests of `crates/kuniumi-cli/tests/cli.rs` compare `metrics.txt` of the
golden world (10,000 cells, seed 1) with these values. Standard library only.
Run from the repository root (see tools/README.md):

    kuniumi generate --seed 1 --cells 10000 --out out/M2-06/10k
    python tools/hypsometry_reference.py out/M2-06/10k
"""

from __future__ import annotations

import json
import math
import struct
import sys
from fractions import Fraction
from pathlib import Path

BIN_M = 100
PERCENTILES = [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]


def read_kfield(path: Path) -> list[float]:
    """The values of a `.kfield` file (kuniumi_io::field_file)."""
    data = path.read_bytes()
    assert data[:8] == b"KUNIFLD\0", path
    (version,) = struct.unpack_from("<I", data, 8)
    assert version == 1, path
    offset = 12
    for _ in range(2):  # name and unit
        (length,) = struct.unpack_from("<I", data, offset)
        offset += 4 + length
    (count,) = struct.unpack_from("<Q", data, offset)
    offset += 8
    return list(struct.unpack_from(f"<{count}d", data, offset))


def main(directory: str) -> None:
    world = Path(directory)
    height = read_kfield(world / "surface_height.kfield")
    area = read_kfield(world / "cell_area.kfield")
    assert len(height) == len(area)
    sea_level = json.loads((world / "world.json").read_text())["sea_level_m"]
    z = Fraction(sea_level)
    elevation = [Fraction(h) - z for h in height]
    weight = [Fraction(a) for a in area]
    total = sum(weight)
    land = [i for i, h in enumerate(height) if not h < sea_level]
    ocean = [i for i, h in enumerate(height) if h < sea_level]

    def mean(cells: list[int], sign: int) -> str:
        if not cells:
            return "none"
        return repr(float(sum(weight[i] * sign * elevation[i] for i in cells) / sum(weight[i] for i in cells)))

    print(f"cells = {len(height)}")
    print(f"sea_level_m = {sea_level!r}")
    print(f"land_fraction = {float(sum(weight[i] for i in land) / total)!r}")
    print(f"mean_land_elevation_m = {mean(land, 1)}")
    print(f"mean_ocean_depth_m = {mean(ocean, -1)}")
    print(f"highest_m = {float(max(elevation))!r}")
    print(f"lowest_m = {float(min(elevation))!r}")

    order = sorted(range(len(height)), key=lambda i: (elevation[i], i))
    running = []
    acc = Fraction(0)
    for i in order:
        acc += weight[i]
        running.append(acc)
    for p in PERCENTILES:
        target = Fraction(p) * total
        k = next(j for j, s in enumerate(running) if s >= target)
        print(f"percentile_{round(p * 100)}_m = {float(elevation[order[k]])!r}")

    bins: dict[int, Fraction] = {}
    for e, w in zip(elevation, weight):
        k = math.floor(e / BIN_M)
        bins[k] = bins.get(k, Fraction(0)) + w
    first, last = min(bins), max(bins)
    share = {k: bins.get(k, Fraction(0)) / total for k in range(first, last + 1)}

    def best(keys: list[int], better) -> int | None:
        found = None
        for k in keys:
            if found is None or better(share[k], share[found]):
                found = k
        return found

    ocean_mode = best([k for k in share if k < 0], lambda a, b: a > b)
    land_mode = best([k for k in share if k >= 0], lambda a, b: a > b)
    valley = None
    if ocean_mode is not None and land_mode is not None:
        valley = best(list(range(ocean_mode + 1, land_mode)), lambda a, b: a < b)
    for name, k in [("ocean_mode", ocean_mode), ("land_mode", land_mode), ("valley", valley)]:
        if k is None:
            print(f"{name}_m = none")
            print(f"{name}_fraction = none")
        else:
            print(f"{name}_m = {(k + 0.5) * BIN_M!r}")
            print(f"{name}_fraction = {float(share[k])!r}")
    separate = valley is not None and share[valley] <= min(share[ocean_mode], share[land_mode]) / 2
    print(f"separate_peaks = {'yes' if separate else 'no'}")
    print(f"bins = {last - first + 1}")
    print(f"first_bin_m = {first * BIN_M}")
    print(f"last_bin_m = {last * BIN_M}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: hypsometry_reference.py WORLD_DIRECTORY")
    main(sys.argv[1])
