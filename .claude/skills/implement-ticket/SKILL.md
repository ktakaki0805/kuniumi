---
name: implement-ticket
description: Implement one kuniumi ticket from docs/tickets end to end — refine it if it is a draft, implement it test-first on a local branch, pass mise run ci, loop with the independent ticket-reviewer subagent until it approves (at most 3 rounds), then squash it into one commit on the local main branch. Nothing is pushed. Use when the user asks to implement, continue or pick up a ticket (an ID such as M1-03, or "the next ticket").
argument-hint: "[ticket-id | next]"
---

# Implement a ticket

Ticket: $ARGUMENTS (an ID such as `M1-03`; empty or `next` means pick the next
one).

Work through the phases in order and track them in the task list. Talk to the
user in Japanese; write code, comments, commit messages and documents in
English (the design document is Japanese).

The rules come from `docs/decisions.md`: tickets and the review loop (D-20),
local integration (D-21), the changelog (D-22), the development environment
(D-23) and dependencies (D-10). In short: the ticket is built and reviewed on a
local branch and squashed into one commit on the local `main`; nothing is ever
pushed, pulled or fetched; every check runs through `mise`.

This skill runs either in a conversation with the owner or inside the
`ticket-implementer` subagent started by `implement-all-tickets`. In the
subagent, "ask the user" and "wait for approval" mean: commit on the ticket
branch and end with a `NEEDS_OWNER` report, as that agent's instructions
describe, and the final report goes to the orchestrator.

## Phase 0 — Pick the ticket

1. With an ID, open `docs/tickets/<ID>-*.md`.
2. Without one, list statuses with `grep -H "^- Status:" docs/tickets/M*.md`,
   walk the order in `docs/tickets/README.md`, and take the first ticket whose
   status is `ready` or `draft` and whose `Depends on:` tickets are all `done`.
   Skip tickets that already have a local `ticket/<ID>-*` branch unless the
   user asked to continue that ticket.
3. If the ticket is `blocked` or `done`, or a dependency is not `done`, stop and
   tell the user why.
4. If a branch for the ticket already exists, switch to it and resume from the
   phase that its commits and the ticket's Result section show.

## Phase 1 — Understand

Read, in this order: `CLAUDE.md`, `docs/conventions.md`, the ticket, the
sections of `docs/design-ja.md` listed in its `Design:` line, the decisions
listed in its `Decisions:` line, the Result sections of its dependency
tickets, and the code it will touch. Write down anything the ticket leaves
undecided.

## Phase 2 — Branch

1. The working tree must be clean; otherwise stop and ask.
2. `mise install` must succeed; otherwise stop and tell the user.
3. `git switch main`, then `git switch -c ticket/<ID>-<slug>`, with the slug
   taken from the ticket's file name. Every later change, including edits to
   the ticket, happens on this branch.

## Phase 3 — Refine

Always for `draft` tickets; for `ready` tickets only where Phase 1 found gaps.

- Turn the Scope and Implementation notes into a concrete plan: public API
  signatures, data structures, the algorithm with formulas and SI units, files
  to touch, tests to write, and how each acceptance criterion will be checked.
- Make every acceptance criterion measurable, and keep the last one:
  "`mise run ci` passes."
- Keep the header lines (`Status`, `Milestone`, `Depends on`, `Crates`,
  `Design`, `Decisions`) accurate; `Design:` entries must be headings of
  `docs/design-ja.md` and `Decisions:` entries must exist in
  `docs/decisions.md`.
- Never settle by assumption anything that needs a decision outside the design
  and `docs/decisions.md`: a borderline or little-known dependency (D-10), a
  parameter the design does not give, a scope change, a dataset (D-16), a
  change to `mise.toml` beyond what the ticket names, or a deviation from a
  cited method. Ask the user, giving the options and your recommendation, and
  stop until they answer.
- For a `draft` ticket, set its status to `ready`, commit the refined ticket
  (`<ID>: refine ticket`), show the user a short summary of the plan, and wait
  for approval before Phase 4.
- If the ticket is too big for one commit (roughly more than 800 changed lines,
  or several independent parts), propose a split and stop. If the user
  approves it, narrow this ticket to the first part, add the rest as new
  `draft` tickets (with dependencies) to `docs/tickets/README.md` in the same
  refine commit, and continue with the narrowed ticket.

## Phase 4 — Implement

- Write tests for the acceptance criteria first where possible, and see them
  fail before the implementation exists.
- Reference values come from an independent source: published values, a hand
  calculation, or a script under `tools/` (Python, D-11), with the script
  referenced in the test comment. Never take expected values from the output
  of the code under test.
- Record golden fingerprints of new outputs only after the invariant tests and
  the visual check pass, and say so in the test comment.
- Stay inside the ticket's scope. Follow-up work becomes a new `draft` ticket
  (copy `docs/tickets/TEMPLATE.md`) listed in `docs/tickets/README.md`.
- Follow `docs/conventions.md`: determinism, SI units with suffixes, the
  dependency rules of D-10, sources cited in doc comments, ported code recorded
  in `THIRD_PARTY.md`.
- Commit small logical steps on the branch, each message starting with
  `<ID>: `. They are squashed later.

## Phase 5 — Self-check

Format, then run the gate, and fix until it passes:

```sh
mise run format
mise run ci      # fmt check, clippy -D warnings, tests (debug, release), MSRV 1.85.0, cargo-deny
```

`mise run check` (without the MSRV build and cargo-deny) is enough while
iterating; the work is done only when `mise run ci` passes. Never fall back to
an unpinned toolchain. Agreement of golden values across platforms cannot be
checked locally; CI confirms it after the owner pushes.

Then:

- For anything spatial, render images into `out/<ID>/` with the CLI or a test
  helper and look at every one with the Read tool. Write down what you see and
  whether it matches the expectation.
- Run the timing measurements the ticket asks for, in release mode.
- For each acceptance criterion, record the evidence: test name, command and
  output excerpt, measured number, or image path with your description. Put it
  in the ticket's Result section as a draft.
- Draft the `CHANGELOG.md` entry under *Unreleased* (Added, Changed, Fixed or
  Removed), plus a *World-breaking* entry if outputs change for existing seeds
  (D-22), so the reviewer can check it.

## Phase 6 — Review loop (at most 3 rounds)

For round r = 1, 2, 3:

1. Commit everything on the branch so the reviewer sees the complete diff
   against `main`.
2. Start a **new** `ticket-reviewer` subagent with the Agent tool (never resume
   an earlier reviewer). Give it facts only:
   - the ticket path and ID, base branch `main`, and the round number r;
   - the acceptance-criteria evidence from the Result section;
   - for r > 1, the previous findings verbatim, each followed by
     "fixed in <commit>" or "not fixed: <reason>".
   Do not argue for the implementation in the prompt.
3. Act on the verdict:
   - `blocker` and `major`: fix them. If you are certain a finding is wrong,
     write the reason, citing the ticket, the design document, the conventions
     or a test.
   - `minor`: fix when cheap; otherwise record it as a follow-up.
   - `nit`: optional.
4. After fixes, run Phase 5 again in full.
5. End the loop on `APPROVE` (no open blocker or major). If blocker or major
   findings are still open after round 3, stop: leave the branch unmerged and
   report to the user with the open findings and your assessment.

## Phase 7 — Land on the local main

1. In the ticket, set `- Status: done` and complete the Result section: what
   was done, deviations from the plan and why, evidence for each criterion,
   measured numbers, a summary of the review rounds, and follow-up tickets.
2. Finalise the `CHANGELOG.md` entry drafted in Phase 5.
3. Update `docs/roadmap.md` if the ticket meets a milestone exit criterion.
   Update `docs/design-ja.md` and `docs/decisions.md` only if the design or a
   decision changed with the user's agreement in Phase 3.
4. Commit these on the branch. If `main` has moved since Phase 2 (for example
   because the owner committed to it), rebase the branch onto `main` and run
   Phase 5 again.
5. Squash into one commit on `main`:
   ```sh
   git switch main
   git merge --squash ticket/<ID>-<slug>
   git commit        # message below
   git branch -D ticket/<ID>-<slug>
   ```
6. Report to the user in Japanese: the ticket, the commit on the local `main`,
   what changed, the key evidence (numbers, images), how many review rounds ran
   and what they caught, follow-up tickets, and anything that needs their
   decision.

### Commit message

```
<ID>: <ticket title>

<one-paragraph summary of what changed and why>

Acceptance criteria: all met; evidence in docs/tickets/<file> (Result).
Review: <n> round(s); <main findings and how they were resolved>.
Follow-ups: <ticket ids, or none>

World-breaking: <no | reason>
```

## Never

- Change a golden value to make a test pass. The only exception is an
  intentional world-breaking change, stated in the commit message and in
  `CHANGELOG.md`.
- Add a dependency outside the D-10 rules, weaken a lint, add `#[allow(...)]`
  to library code, or mark a failing test `#[ignore]`. Ask the user instead.
- Touch files unrelated to the ticket, apart from new follow-up tickets and the
  documents named above.
- Skip the review loop, land work it has not approved, or report a check as
  passing without having run it.
- Run `git push`, `git pull` or `git fetch`, add or change remotes, or rewrite
  the history of `main`.
