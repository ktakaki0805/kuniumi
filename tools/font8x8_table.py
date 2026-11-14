"""Glyph tables of the bitmap font of `kuniumi-viewer` and the hash its test quotes.

Ticket MV-06. Uses only the standard library. Downloads three files of
font8x8 by Daniel Hepper (after Marcel Sondaar's IBM VGA fonts; public
domain), <https://github.com/dhepper/font8x8>, at a pinned commit and checks
their SHA-256:

- `font8x8_basic.h`, the glyphs of U+0000 to U+007F;
- `font8x8_ext_latin.h`, the glyphs of U+00A0 to U+00FF;
- `README`, whose "Encoding" section documents the glyph of `A`.

Prints:

1. `BASIC` (U+0020 to U+007E) and `DEGREE` (U+00B0) as the Rust source of
   `crates/kuniumi-viewer/src/font.rs`, with the source's comment on each
   glyph, for comparison with that file;
2. the FNV-1a 64-bit hash of the 768 bytes of `BASIC` followed by `DEGREE`,
   quoted by the test `the_table_matches_its_source` in `font.rs`.

It also checks that the README's `A` is the table's.

Usage, from the repository root:

    mise exec -- python tools/font8x8_table.py
"""

import hashlib
import re
import urllib.request

BASE_URL = (
    "https://raw.githubusercontent.com/dhepper/font8x8/"
    "8e279d2d864e79128e96188a6b9526cfa3fbfef9/"
)
FILES = {
    "font8x8_basic.h": "49d8df366296b203ca3211bc0672cf2a762135bf12710735b6292756b19dffd5",
    "font8x8_ext_latin.h": "e0582e9a6b3f85dbb983c72d95c9a50c3547ad5e895d16ef167e20777d053f1b",
    "README": "885cc022fa9c844695b364e77e6a18ab2be74678a181827f1ce893071205cb65",
}

GLYPH = re.compile(r"\{((?:\s*0x[0-9A-Fa-f]{2}\s*,){7}\s*0x[0-9A-Fa-f]{2})\s*\},?\s*(// .*)$")


def download(name):
    with urllib.request.urlopen(BASE_URL + name) as response:
        data = response.read()
    digest = hashlib.sha256(data).hexdigest()
    if digest != FILES[name]:
        raise SystemExit(f"{name}: SHA-256 {digest}, expected {FILES[name]}")
    return data.decode("utf-8")


def glyphs(src, first, count):
    """The glyphs of a header's table, by position: code point `first + i`
    for entry `i`, as the C array is indexed -> (8 bytes, the source's
    comment). The comments are kept for reading only: some of
    `font8x8_ext_latin.h` name the wrong code point (U+00B3 and U+00B4 are
    labelled U+00B2)."""
    entries = []
    for line in src.splitlines():
        m = GLYPH.search(line)
        if m:
            rows = [int(v, 16) for v in m.group(1).split(",")]
            entries.append((rows, m.group(2).rstrip()))
    assert len(entries) == count, f"{len(entries)} glyphs, expected {count}"
    return {first + i: entry for i, entry in enumerate(entries)}


def fnv1a64(data):
    h = 0xCBF29CE484222325
    for b in data:
        h ^= b
        h = (h * 0x100000001B3) % (1 << 64)
    return h


def rust_row(rows, comment):
    return "    [" + ", ".join(f"0x{r:02X}" for r in rows) + "], " + comment


def main():
    basic = glyphs(download("font8x8_basic.h"), 0x00, 0x80)
    latin = glyphs(download("font8x8_ext_latin.h"), 0xA0, 0x60)
    readme = download("README")
    for code in range(0x20, 0x7F):
        assert f"// U+{code:04X} " in basic[code][1] + " ", basic[code]
    assert latin[0xB0][1] == "// U+00B0 (degree)", latin[0xB0]

    documented = re.search(r"is encoded as\s*\n\s*\{([^}]*)\}", readme).group(1)
    a = [int(v, 16) for v in documented.split(",")]
    assert a == basic[0x41][0], f"README's A {a} differs from the table's {basic[0x41][0]}"

    print("const BASIC: [[u8; 8]; 95] = [")
    data = bytearray()
    for code in range(0x20, 0x7F):
        rows, comment = basic[code]
        print(rust_row(rows, comment))
        data.extend(rows)
    print("];")
    rows, comment = latin[0xB0]
    print(f"const DEGREE: [u8; 8] = [{', '.join(f'0x{r:02X}' for r in rows)}]; {comment}")
    data.extend(rows)
    assert len(data) == 768
    print(f"// FNV-1a 64 of BASIC then DEGREE ({len(data)} bytes): 0x{fnv1a64(data):016x}")


if __name__ == "__main__":
    main()
