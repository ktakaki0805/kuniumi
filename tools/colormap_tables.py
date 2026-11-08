"""Colour tables of `kuniumi_io::colormap` and reference values for its tests.

Ticket M1-08. Uses only the standard library. Downloads two files at pinned
commits and checks their SHA-256:

- `colormaps.py` of the BIDS colormap project (viridis by Nathaniel J. Smith,
  Stefan van der Walt and Eric Firing; CC0), whose `_viridis_data` is the
  256-entry viridis table that matplotlib ships;
- `colorbrewer_schemes.js` of ColorBrewer (Cynthia Brewer, Mark Harrower and
  The Pennsylvania State University; Apache-2.0).

Prints:

1. `VIRIDIS`: every 5th of the 256 viridis entries (52 control points at
   t = k/51), as Rust source, and the largest deviation of linear
   interpolation through them from linear interpolation through all 256
   entries, in 8-bit units;
2. reference viridis colours at t = j/20, interpolated linearly through all
   256 entries and rounded half up to 8 bits (the test
   `viridis_matches_the_full_table` in `colormap/tests.rs`);
3. `DIVERGING` (RdBu, 11 classes, reversed so that blue is low) and
   `CATEGORICAL` (Paired, 12 classes), as Rust source;
4. the smallest CIE76 colour difference between two colours of Paired
   (sRGB, D65 white), quoted by the test `categorical_colours_are_distinct`.

Usage, from the repository root:

    mise exec -- python tools/colormap_tables.py
"""

import ast
import hashlib
import re
import urllib.request

BIDS_URL = (
    "https://raw.githubusercontent.com/BIDS/colormap/"
    "bc549477db0c12b54a5928087552ad2cf274980f/colormaps.py"
)
BIDS_SHA256 = "b7478d10e2ccfdb22128867e137f16f55b687e9c225065e91bcbd6bb9e6c0352"
BREWER_URL = (
    "https://raw.githubusercontent.com/axismaps/colorbrewer/"
    "7d135fc4e19eda73f2eb1bf55fcdf4a04fe4881f/colorbrewer_schemes.js"
)
BREWER_SHA256 = "52e6a6edbbe4bbb95467cb84fb7db8fcbb1d6d219c975851c4c9d79334f00a41"

STEP = 5  # every 5th entry: 0, 5, ..., 255


def download(url, sha256):
    with urllib.request.urlopen(url) as response:
        data = response.read()
    digest = hashlib.sha256(data).hexdigest()
    if digest != sha256:
        raise SystemExit(f"{url}: SHA-256 {digest}, expected {sha256}")
    return data.decode("utf-8")


def viridis_table(src):
    start = src.index("_viridis_data = ") + len("_viridis_data = ")
    end = src.index("]]", start) + 2
    table = ast.literal_eval(src[start:end])
    assert len(table) == 256
    return table


def brewer(src, name, classes):
    block = src[src.index(name + ":") :]
    block = block[: block.index("}")]
    found = re.search(r"\b%d: \[([^\]]*)\]" % classes, block).group(1)
    return [tuple(int(v) for v in m) for m in re.findall(r"rgb\((\d+),(\d+),(\d+)\)", found)]


def lerp(points, t):
    """Linear interpolation through equally spaced points, t in [0, 1]."""
    s = t * (len(points) - 1)
    k = min(int(s), len(points) - 2)
    f = s - k
    return [a + f * (b - a) for a, b in zip(points[k], points[k + 1])]


def round_half_up(x):
    return int(x + 0.5)


def lab(rgb):
    """CIELAB (D65) of an 8-bit sRGB colour (IEC 61966-2-1, CIE 15)."""

    def linear(c):
        c /= 255.0
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4

    r, g, b = (linear(c) for c in rgb)
    x = 0.4124 * r + 0.3576 * g + 0.1805 * b
    y = 0.2126 * r + 0.7152 * g + 0.0722 * b
    z = 0.0193 * r + 0.1192 * g + 0.9505 * b

    def f(t):
        return t ** (1 / 3) if t > (6 / 29) ** 3 else t / (3 * (6 / 29) ** 2) + 4 / 29

    fx, fy, fz = f(x / 0.95047), f(y), f(z / 1.08883)
    return (116 * fy - 16, 500 * (fx - fy), 200 * (fy - fz))


def main():
    table = viridis_table(download(BIDS_URL, BIDS_SHA256))
    points = table[::STEP]
    assert len(points) == 52 and points[-1] == table[255]

    print(f"// VIRIDIS: {len(points)} control points")
    for p in points:
        print(f"    [{p[0]:.6f}, {p[1]:.6f}, {p[2]:.6f}],")
    worst = 0.0
    for i in range(255 * 20 + 1):
        t = i / (255 * 20)
        a, b = lerp(points, t), lerp(table, t)
        worst = max(worst, max(abs(x - y) for x, y in zip(a, b)) * 255)
    print(f"// largest deviation from the 256-entry table: {worst:.3f} (8-bit units)")

    print("// viridis reference colours at t = j/20")
    for j in range(21):
        c = [round_half_up(255 * v) for v in lerp(table, j / 20)]
        print(f"    ({j}, [{c[0]}, {c[1]}, {c[2]}]),")

    src = download(BREWER_URL, BREWER_SHA256)
    print("// DIVERGING: ColorBrewer RdBu, 11 classes, reversed")
    for c in reversed(brewer(src, "RdBu", 11)):
        print(f"    [{c[0]}, {c[1]}, {c[2]}],")
    paired = brewer(src, "Paired", 12)
    print("// CATEGORICAL: ColorBrewer Paired, 12 classes")
    for c in paired:
        print(f"    [{c[0]}, {c[1]}, {c[2]}],")
    labs = [lab(c) for c in paired]
    closest = min(
        (sum((u - v) ** 2 for u, v in zip(labs[i], labs[j])) ** 0.5, i, j)
        for i in range(12)
        for j in range(i + 1, 12)
    )
    print(f"// smallest CIE76 difference in Paired: {closest[0]:.2f} (classes {closest[1]}, {closest[2]})")


if __name__ == "__main__":
    main()
