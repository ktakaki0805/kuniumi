# M2-09 Rustdoc warnings

- Status: done
- Milestone: M2
- Depends on: none
- Crates: kuniumi-io, kuniumi-sphere
- Design: プロジェクト方針
- Decisions: D-23, D-25

## Goal

`cargo doc --workspace --no-deps` builds without warnings, and the owner has
decided whether `mise run ci` keeps it that way.

## Scope

In:
- The warning the ticket was opened for: the public documentation of
  `kuniumi_io::colormap` links to the private item `VIRIDIS`
  (`crates/kuniumi-io/src/colormap.rs`, line 7; lint
  `rustdoc::private_intra_doc_links`), so
  `RUSTDOCFLAGS="-D warnings" cargo doc -p kuniumi-io` fails.
- A second warning that the verification command finds with the pinned
  toolchain (Rust 1.97.0): the module documentation of
  `kuniumi_sphere::points` writes `` [`TAU`](core::f64::consts::TAU) ``, and
  `TAU` is already in scope through `use core::f64::consts::{PI, TAU}`
  (`crates/kuniumi-sphere/src/points.rs`, line 36; lint
  `rustdoc::redundant_explicit_links`). It dates from M1-04. The Goal covers
  the whole workspace, so it is fixed here too.
- With the owner: whether a `doc` task (`cargo doc --workspace --no-deps`
  with `-D warnings`) joins `mise run check` (a change to `mise.toml`, D-23)
  and the CI workflow (D-25).

Out:
- Rewriting documentation beyond the warnings.
- Making `VIRIDIS` public: it is an implementation detail of
  `Colormap::Viridis`.
- Documenting private items (`--document-private-items`): the check is on
  the public documentation, which is what users of the crates read.
- `CLAUDE.md` and the files under `.claude/`: they describe the checks in
  passing (`mise run check  # fmt check, clippy -D warnings, tests …`); they
  are instruction files that the owner edits, so if a `doc` task is added
  their wording is listed in the Result for the owner instead of changed
  here.

## Implementation notes

Found by the review of M2-02 (R2-3); the link predates M2-02, which does not
touch the file.

### The fixes (documentation only)

1. `crates/kuniumi-io/src/colormap.rs`, line 7: `` [`VIRIDIS`] `` becomes the
   code span `` `VIRIDIS` `` (no link). The sentence stays as it is; the
   constant keeps its own doc comment, which points back to the module
   documentation.
2. `crates/kuniumi-sphere/src/points.rs`, line 36:
   `` [`TAU`](core::f64::consts::TAU) `` becomes `` [`TAU`] ``, which rustdoc
   resolves through the `use` to the same item, `core::f64::consts::TAU`.

Neither change touches code, so no output changes (World-breaking: no) and
every golden test passes unchanged.

### The `doc` task (only if the owner says yes)

`mise.toml`:

```toml
[tasks.doc]
description = "Build the documentation, warnings are errors"
env = { RUSTDOCFLAGS = "-D warnings" }
run = "cargo doc --workspace --no-deps --locked"

[tasks.check]
description = "fmt, clippy, doc and tests in debug and release"
run = [{ task = "fmt" }, { task = "clippy" }, { task = "doc" }, { task = "test" }]
```

It runs with the development toolchain only (not the MSRV one); new rustdoc
lints arrive when the toolchain is raised on purpose, as clippy lints do
(D-23). `.github/workflows/ci.yml`: the `lint` job runs `mise run doc` after
`mise run clippy`, so CI keeps running the same tasks as `mise run ci`
(D-25). `docs/decisions.md`: D-23 says that `mise run check` also builds the
documentation with warnings as errors, and D-25 adds "rustdoc" to the Linux
checks. `docs/conventions.md` (Code style) and the commands of `README.md`:
the list of what `mise run ci` runs gains the documentation build.
`CHANGELOG.md`: an *Added* entry for the task and a *Fixed* entry for the
two links.

Cost, measured on this machine (24 threads): about 3.1 s from an empty target
directory, under 1 s when the dependencies are already checked.

If the owner says no, no file besides the two sources, this ticket and
`CHANGELOG.md` changes, and the decision is recorded in the Result.

### Tests

There is no unit test to write for a documentation lint; the test is the
documentation build itself. Before the fixes, the verification command must
fail on exactly the two warnings above; after them it must pass. If the
`doc` task is added, it is first run against the unfixed sources to see it
fail.

## Acceptance criteria

- [x] Before the fixes, `RUSTDOCFLAGS="-D warnings" mise exec -- cargo doc
      --workspace --no-deps --locked` fails, and without `-D warnings` the
      build reports exactly two warnings: `colormap.rs:7`
      (`private_intra_doc_links`) and `points.rs:36`
      (`redundant_explicit_links`). After them, the same command exits 0 and
      prints no `warning`.
- [x] In the generated HTML, the `TAU` link of `kuniumi_sphere::points`
      still points to `core::f64::consts::TAU` (the `href` in
      `target/doc/kuniumi_sphere/points/index.html`), and `VIRIDIS` in
      `target/doc/kuniumi_io/colormap/index.html` is code text without a
      link.
- [x] The owner's decision on a `doc` task is recorded in the Result. If
      yes: `mise run doc` exists, fails on the unfixed sources and passes
      after the fixes; `mise run check` (and so `mise run ci`) runs it; the
      CI `lint` job runs it; D-23, D-25 and the conventions say so.
- [x] No output changes: the golden tests pass unchanged.
- [x] `mise run ci` passes.

## Verification

```sh
RUSTDOCFLAGS="-D warnings" mise exec -- cargo doc --workspace --no-deps --locked
grep -o 'href="[^"]*TAU[^"]*"' target/doc/kuniumi_sphere/points/index.html
grep -c 'VIRIDIS' target/doc/kuniumi_io/colormap/index.html
mise run ci
```

## Result

Owner's decision (2026-11-12, asked through `implement-all-tickets`):
option A, "doc を check/ci に含める": a `doc` task in `mise.toml`, run by
`mise run check` (and so `mise run ci`) and by the CI `lint` job, with
D-23, D-25 and the conventions updated.

What was done:
- `crates/kuniumi-io/src/colormap.rs`, line 7: `` [`VIRIDIS`] `` is now the
  code span `` `VIRIDIS` ``.
- `crates/kuniumi-sphere/src/points.rs`, line 36:
  `` [`TAU`](core::f64::consts::TAU) `` is now `` [`TAU`] ``; the word
  "before" moved up from line 37 so the paragraph keeps its width.
- `mise.toml`: the task `doc` (`cargo doc --workspace --no-deps --locked`,
  `RUSTDOCFLAGS = "-D warnings"`), run by `check` between `clippy` and
  `test`.
- `.github/workflows/ci.yml`: the `lint` job runs `mise run doc` after
  `mise run clippy`.
- `docs/decisions.md` (D-23, D-25), `docs/conventions.md` (Code style),
  `README.md` (commands) and `CHANGELOG.md` (*Added*: the task; *Fixed*:
  the two links).

Not changed, for the owner: three instruction files still list the checks
without the documentation build. `CLAUDE.md` line 26
(`mise run check  # fmt check, clippy -D warnings, tests in debug and
release`), `.claude/skills/implement-ticket/SKILL.md` line 109 (the comment
after `mise run ci`) and `.claude/agents/ticket-reviewer.md` line 30 ("Run
`mise run ci` yourself (fmt, clippy, tests in debug and release, MSRV,
cargo-deny …)"). Each would gain "doc" (or "rustdoc `-D warnings`").

Evidence:
1. Before the fixes (commit "doc task in mise check, failing on the two
   warnings"), `mise run doc` exited 101 with
   `error: redundant explicit link target --> crates/kuniumi-sphere/src/points.rs:36:54`.
   Whether the second error,
   `error: public documentation for `colormap` links to private item `VIRIDIS` --> crates/kuniumi-io/src/colormap.rs:7:20`,
   is printed too depends on job scheduling, since cargo stops after the
   first failing crate (it was in this run, not in the reviewer's);
   `-p kuniumi-io` with `-D warnings` always fails on it. The same build
   without `-D warnings` reported exactly these two warnings (one per
   crate). After the fixes,
   `RUSTDOCFLAGS="-D warnings" mise exec -- cargo doc --workspace --no-deps
   --locked` exits 0, and the same build without `-D warnings` prints 0
   lines containing `warning`.
2. `target/doc/kuniumi_sphere/points/index.html` has
   `<a href="https://doc.rust-lang.org/1.97.0/core/f64/consts/constant.TAU.html" title="constant core::f64::consts::TAU">`
   both before the fix (built from the stashed sources) and after it.
   `target/doc/kuniumi_io/colormap/index.html` has
   `under CC0. <code>VIRIDIS</code> holds every 5th of the`, and 0
   occurrences of `constant.VIRIDIS`.
3. Decision recorded above. `mise tasks` lists
   `doc  Build the documentation, warnings are errors` and
   `check  fmt, clippy, doc and tests in debug and release`; the
   `mise run ci` log shows `[doc] $ cargo doc --workspace --no-deps
   --locked` between `[clippy]` and `[test-debug]`, documenting all seven
   crates (kuniumi-math, -sphere, -world, -sim, -io, -config, -cli).
   The `lint` job of `ci.yml` has the step `- run: mise run doc`.
4. The diff changes no code, only two doc comments; no test file and no
   golden value changes, and every golden test passes in debug and release
   within `mise run ci`.
5. `mise run ci` exits 0 (1 min 40 s): fmt, clippy `-D warnings`, doc
   `-D warnings`, 921 tests over debug and release (0 failed; 4 ignored per
   profile, the existing timing tests), MSRV 1.85.0, cargo-deny (advisories,
   bans, licences, sources ok; its warnings about unused licence allowances
   and the duplicate `miniz_oxide` predate this ticket).

Cost of the `doc` task on this machine (24 threads): 3.1 s from an empty
target directory, 0.4–0.6 s when the dependencies are already checked.

Deviations from the draft: the draft named one warning; the verification
command on the whole workspace found a second (`points.rs:36`, from M1-04),
fixed here because the Goal covers the workspace, so `Crates:` gained
`kuniumi-sphere`. `README.md` and the toolchain comment in `mise.toml` were
updated besides the files the plan named, because they describe the same
checks.

Review:
- Round 1 (APPROVE): the reviewer re-ran `mise run ci`, rebuilt the
  documentation of commit 6d046e5 from a `git archive` tree (the task fails;
  exactly two warnings without `-D warnings`) and of HEAD in fresh target
  directories, and checked both HTML pages. Two nits: R1-1, the comment on
  the toolchain pin in `mise.toml` did not mention rustdoc lints as D-23 now
  does; R1-2, the evidence said "4 ignored" (4 per profile, 8 in all) and
  that one `mise run doc` run printed both errors, which depends on job
  scheduling. Both fixed.
- Round 2 (APPROVE, no findings): the reviewer re-ran `mise run ci` and the
  builds of 6d046e5 and HEAD, and confirmed both fixes.

Follow-ups: none. The three instruction files listed above are for the
owner to update.
