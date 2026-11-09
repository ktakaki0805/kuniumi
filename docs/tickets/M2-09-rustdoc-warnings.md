# M2-09 Rustdoc warnings

- Status: draft
- Milestone: M2
- Depends on: none
- Crates: kuniumi-io
- Design: プロジェクト方針
- Decisions: D-23

## Goal

`cargo doc --workspace --no-deps` builds without warnings, and the owner has
decided whether `mise run ci` keeps it that way.

## Scope

In:
- The existing warning: the public documentation of
  `kuniumi_io::colormap` links to the private item `VIRIDIS`
  (`crates/kuniumi-io/src/colormap.rs`, line 7), so
  `RUSTDOCFLAGS="-D warnings" cargo doc -p kuniumi-io` fails.
- With the owner: whether a `doc` task (`cargo doc --workspace --no-deps`
  with `-D warnings`) joins `mise run check` (a change to `mise.toml`, D-23).

Out:
- Rewriting documentation beyond the warnings.

## Implementation notes

Found by the review of M2-02 (R2-3); the link predates M2-02, which does not
touch the file.

## Acceptance criteria

- [ ] `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps` passes.
- [ ] The owner's decision on a `doc` task is recorded, and applied if yes.
- [ ] `mise run ci` passes.

## Verification

```sh
RUSTDOCFLAGS="-D warnings" mise exec -- cargo doc --workspace --no-deps
```

## Result

_(filled in when the ticket is implemented)_
