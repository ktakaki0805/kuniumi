# M2-16 Crate lists of the decisions and the design

- Status: done
- Milestone: M2
- Depends on: M2-07
- Crates: none (documentation)
- Design: プロジェクト方針
- Decisions: D-06

## Goal

The crate lists of `docs/decisions.md` (D-06) and `docs/design-ja.md`
(プロジェクト方針, row 言語・構成) name `kuniumi-config`, which exists since
M2-07, as an existing crate, as `CLAUDE.md` and `README.md` already do.

## Scope

In:
- D-06: a `kuniumi-config` row in the crate table, between `kuniumi-world`
  and `kuniumi-sim` (the order of `CLAUDE.md`, `README.md` and the
  dependencies: it depends on `kuniumi-math` and `kuniumi-world`, and only
  `kuniumi-cli` depends on it), and `kuniumi-config` out of the sentence on
  the planned crates.
- With the owner: whether D-06 also gets the note
  "(changed YYYY-MM-DD: what changed)" that the head of `docs/decisions.md`
  asks for when a decision changes. M2-08 (D-12) and M2-09 (D-23, D-25)
  edited decisions in place without that note, naming their ticket in the
  text instead ("chosen with the owner in M2-08", "added by M2-09").
- The design's プロジェクト方針 table: the 言語・構成 row reads
  `Rust の cargo workspace（math / sphere / world / config / sim / io / cli）`.
- `CHANGELOG.md`: a *Fixed* entry (D-22: every ticket's commit updates it).

Out:
- Any change to the decisions themselves: the rule that dependencies point
  downwards, the sentence that `kuniumi-sim` never depends on `kuniumi-io`,
  and the other planned crates stay as they are.
- The date in the head of `docs/decisions.md` ("as of 2026-11-08"): M2-08
  and M2-09 changed decisions after it and left it; changing what it means
  is not this ticket's question.
- `CLAUDE.md` and `README.md`: they already list `kuniumi-config` in this
  order.
- The design's 技術スタック section: it already names toml for the settings
  files and points to D-10 for the list of dependencies.

## Implementation notes

Found by the review of M2-07 (R1-1). The `implement-ticket` skill changes
`docs/decisions.md` and `docs/design-ja.md` only with the owner's agreement,
so M2-07 left them as they were.

### D-06 after the change

```markdown
| Crate | Role |
| --- | --- |
| `kuniumi-math` | deterministic math, random numbers, vectors, SI units |
| `kuniumi-sphere` | the spherical cell grid |
| `kuniumi-world` | planet parameters, settings, the world data model |
| `kuniumi-config` | settings files (TOML, D-15) |
| `kuniumi-sim` | the generation pipeline (no I/O) |
| `kuniumi-io` | file formats, images, fingerprints |
| `kuniumi-cli` | the `kuniumi` command |

`kuniumi-sim` never depends on `kuniumi-io`. Planned crates: `kuniumi-viewer`
(D-19), `kuniumi-gpu` (D-08), `kuniumi-ffi` (engine bindings, D-07).
```

If the owner wants the change note, the paragraph ends with
"(changed 2026-11-13: `kuniumi-config`, built by M2-07, moved from the
planned crates into the table)", the date being the day the ticket lands.
Otherwise the text stays as above, as the in-place edits of M2-08 and M2-09
did.

### The design after the change

```markdown
| 言語・構成 | Rust の cargo workspace（math / sphere / world / config / sim / io / cli） | D-06 |
```

### Changelog

Under *Fixed*: "D-06 (`docs/decisions.md`) and the project table of the
design list `kuniumi-config` as an existing crate, no longer as planned
(M2-16)."

### Tests

There is no code to test. The criteria are checked with `grep` and a shell
loop that compares the D-06 table with the directories under `crates/`; the
loop is first run against `main` to see it report `kuniumi-config` missing.
No output changes (World-breaking: no); the golden tests pass unchanged.

## Acceptance criteria

- [x] The D-06 table has exactly the seven rows `kuniumi-math`,
      `kuniumi-sphere`, `kuniumi-world`, `kuniumi-config`, `kuniumi-sim`,
      `kuniumi-io`, `kuniumi-cli`, in this order, which are the crates under
      `crates/` (each `crates/<name>/Cargo.toml` has `name = "<name>"`) in
      the order of `CLAUDE.md`'s workspace list; on `main` the same check
      reports `kuniumi-config` missing.
- [x] The planned crates of D-06 are `kuniumi-viewer`, `kuniumi-gpu` and
      `kuniumi-ffi`, with their decisions as before; `kuniumi-config` is not
      among them.
- [x] D-06 ends with the change note if and only if the owner chose it; the
      owner's choice is recorded in the Result.
- [x] The 言語・構成 row of the design's プロジェクト方針 table reads
      `math / sphere / world / config / sim / io / cli`.
- [x] No other line of `docs/decisions.md` or `docs/design-ja.md` changes
      (`git diff main --stat` and the diff itself).
- [x] `CHANGELOG.md` has the *Fixed* entry.
- [x] `mise run ci` passes.

## Verification

```sh
grep -n "kuniumi-config" docs/decisions.md docs/design-ja.md
grep -n "言語・構成" docs/design-ja.md
sed -n '/^\*\*D-06/,/^\*\*D-07/p' docs/decisions.md
# the table against the crates
diff <(sed -n '/^\*\*D-06/,/^\*\*D-07/p' docs/decisions.md \
         | grep -o '^| `kuniumi-[a-z]*`' | tr -d '|` ') \
     <(sed -n '/^## Workspace/,/^## /p' CLAUDE.md | grep -o '^kuniumi-[a-z]*')
for d in crates/*/; do
  n=$(basename "$d")
  grep -q "^name = \"$n\"" "$d/Cargo.toml" \
    && sed -n '/^\*\*D-06/,/^\*\*D-07/p' docs/decisions.md | grep -q "^| \`$n\`" \
    || echo "missing: $n"
done
git diff main --stat
mise run ci
```

## Result

Owner's decision (asked after the refine commit): approve the plan and end
D-06 with the change note that the head of `docs/decisions.md` asks for.

What was done:
- `docs/decisions.md`, D-06: the row `` | `kuniumi-config` | settings files
  (TOML, D-15) | `` between `kuniumi-world` and `kuniumi-sim`; the planned
  crates are `kuniumi-viewer` (D-19), `kuniumi-gpu` (D-08) and
  `kuniumi-ffi` (engine bindings, D-07); the paragraph ends with
  "(changed 2026-11-13: `kuniumi-config`, built by M2-07, moved from the
  planned crates into the table)".
- `docs/design-ja.md`, プロジェクト方針: the 言語・構成 row reads
  `Rust の cargo workspace（math / sphere / world / config / sim / io / cli）`.
- `CHANGELOG.md`: the *Fixed* entry of the plan.

Deviations from the plan: none.

Evidence:

1. The table against the crates: on `main`, the `diff` of the D-06 table
   with `CLAUDE.md`'s workspace list prints
   `3a4 > kuniumi-config` (exit 1) and the loop over `crates/*/` prints
   `missing: kuniumi-config`. After the change the `diff` prints nothing
   (exit 0), so the table has the seven rows in `CLAUDE.md`'s order, and the
   loop prints nothing: every directory under `crates/` has a
   `Cargo.toml` with its own `name` and a row in the table. The order
   follows the dependencies: `kuniumi-config` depends on `kuniumi-math` and
   `kuniumi-world` (its `Cargo.toml`), and only `kuniumi-cli` depends on it.
2. Planned crates: `sed -n '/^\*\*D-06/,/^\*\*D-07/p' docs/decisions.md`
   shows "Planned crates: `kuniumi-viewer` (D-19), `kuniumi-gpu` (D-08),
   `kuniumi-ffi` (engine bindings, D-07)."; `grep -n kuniumi-config
   docs/decisions.md` gives lines 40 (the table row), 47 (the change note),
   94 (the D-10 table) and 135 (D-15), none in the planned crates.
3. Change note: D-06 ends with the note above, as the owner chose.
4. Design: `grep -n "言語・構成" docs/design-ja.md` gives line 26 with
   `math / sphere / world / config / sim / io / cli`.
5. No other line changes: `git diff main --stat` lists `CHANGELOG.md` (+2),
   `docs/decisions.md` (8 lines: the new row, the planned-crates sentence
   and the note), `docs/design-ja.md` (1 line) and this ticket.
6. `CHANGELOG.md`: under *Fixed*, "D-06 (`docs/decisions.md`) and the
   project table of the design list `kuniumi-config` as an existing crate,
   no longer as planned (M2-16)."
7. `mise run ci` passes: fmt, clippy `-D warnings`, doc, 496 tests in debug
   and 487 in release (4 ignored in each, as before), MSRV 1.85.0,
   cargo-deny (advisories, bans, licences, sources ok; warnings only for
   allowed licences not in use and the known duplicate `miniz_oxide`).

No image: the ticket changes documentation only, nothing spatial. No code
changes, so no output changes (World-breaking: no).

Review:
- Round 1 (APPROVE, no findings): the reviewer re-ran `mise run ci` (496
  and 487 tests, as above), the Verification checks on the branch and
  against `main`, checked the dependencies of `kuniumi-config` in the
  `Cargo.toml` files, searched `docs/`, `README.md` and `CLAUDE.md` for
  other stale crate lists (none), and confirmed that the change note
  follows the format at the head of `docs/decisions.md` and that only D-06
  and row 26 of the design change.

Follow-ups: none.
