# M2-16 Crate lists of the decisions and the design

- Status: draft
- Milestone: M2
- Depends on: M2-07
- Crates: none (documentation)
- Design: プロジェクト方針
- Decisions: D-06

## Goal

The crate lists of `docs/decisions.md` (D-06) and `docs/design-ja.md`
(プロジェクト方針, 言語・構成) name `kuniumi-config`, which exists since
M2-07, as `CLAUDE.md` and `README.md` already do.

## Scope

In:
- D-06: `kuniumi-config` (settings files, D-15) in the crate table, out of
  the planned crates, with the "(changed YYYY-MM-DD: …)" note the decisions
  file asks for, if the owner agrees that this is the right way to record it.
- The design's プロジェクト方針 table: the workspace row reads
  math / sphere / world / config / sim / io / cli.

Out:
- Any change to the decisions themselves.

## Implementation notes

Found by the review of M2-07 (R1-1). The `implement-ticket` skill changes
`docs/decisions.md` and `docs/design-ja.md` only with the owner's agreement,
so M2-07 left them as they were.

## Acceptance criteria

- [ ] D-06 and the design's table list `kuniumi-config` as an existing crate.
- [ ] `mise run ci` passes.

## Verification

```sh
grep -n "kuniumi-config" docs/decisions.md docs/design-ja.md
```

## Result

_(filled in when the ticket is implemented)_
