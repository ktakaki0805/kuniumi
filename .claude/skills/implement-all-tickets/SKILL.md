---
name: implement-all-tickets
description: Implement kuniumi tickets one after another in the order of docs/tickets/README.md, each in a fresh ticket-implementer subagent that follows the implement-ticket skill (review loop, one squashed commit on the local main, nothing pushed), until a target is reached or the owner is needed. Started only by the owner, to implement all, the remaining or several tickets, or a whole milestone.
argument-hint: "[<milestone> | <ticket-id> | <count>]"
disable-model-invocation: true
---

# Implement tickets in sequence

Target: $ARGUMENTS

You are the orchestrator. You do not implement tickets yourself: each ticket is
implemented by a fresh `ticket-implementer` subagent, which follows the
`implement-ticket` skill completely (refinement, tests first, `mise run ci`,
the `ticket-reviewer` loop, one squashed commit on the local `main`). Your job
is to choose tickets, relay the owner's decisions, verify each result, and stop
at the right moment. Talk to the user in Japanese.

Every rule of `implement-ticket` and `docs/decisions.md` still holds for every
ticket — tickets and the review loop (D-20), local integration (D-21), the
changelog (D-22) and the development environment (D-23). Running many tickets
in a row never relaxes them. Nothing is pushed, pulled or fetched.

## Phase 0 — Target

Work out the target set from the argument:

| Argument | Target set | The run ends when |
| --- | --- | --- |
| empty or `all` | every ticket | no ticket is left to start |
| a milestone, e.g. `M3` or `MV` | that milestone's tickets and everything they depend on, transitively | every ticket of the milestone is `done` |
| a ticket ID, e.g. `M2-05` | that ticket and everything it depends on, transitively | that ticket is `done` |
| a number N | every ticket | N tickets have landed in this run |

Dependencies come from the `Depends on:` lines; statuses only from each
ticket's `- Status:` line. If the argument matches none of these, ask the user
what they meant.

## Phase 1 — Preflight

1. Read `CLAUDE.md` and `docs/tickets/README.md`.
2. Git: the working tree must be clean and `main` checked out; otherwise stop
   and ask.
3. Interrupted work: if a local `ticket/<ID>-*` branch exists and that ticket
   is in the target set, it is resumed first (its implementer continues from
   the branch). If it is outside the target set, stop and ask.
4. Baseline: `mise install`, then `mise run ci` on `main`. If either fails,
   stop and report; never build on a broken `main` or an unpinned toolchain.
5. Plan: list, in README order, the tickets of the target set that are not
   `done`, and how many of them are `draft` (each draft stops once for the
   owner to approve its refined plan). Send the plan with SendUserMessage and
   create one task per ticket in the task list.
6. Run log: create `out/implement-all/<YYYYMMDD-HHMM>.md` (`out/` is
   git-ignored) with the argument, the plan and the starting commit of `main`.
   Append to it after every event below. If this conversation is compacted,
   re-read the run log and `git log --oneline` before continuing.

## Phase 2 — Loop

Repeat until the run ends:

1. **Pick.** The next ticket is the first one in README order that is in the
   target set, not `done`, and whose dependencies are all `done`. If none is
   left but the target is not reached (only `blocked` tickets, or tickets
   waiting on them), the run ends; report what blocks it.
2. **Delegate.** Mark the ticket's task in progress and start a **new**
   `ticket-implementer` subagent in the foreground (`run_in_background: false`)
   with only: the ticket ID and file, "implement this ticket by following the
   implement-ticket skill", and, if it was interrupted earlier, that its branch
   exists. Never reuse an implementer for a second ticket.
3. **Handle the report** by its `OUTCOME`:
   - `NEEDS_OWNER`: ask the user with AskUserQuestion, using the report's
     question and options (recommended option first; the question names the
     ticket and stands on its own). For a draft plan, put the plan summary in
     the question and offer "承認して続行" and "ここで中止"; requested changes
     arrive as the user's free-text answer. Send the answer verbatim to the
     **same** implementer with SendMessage and wait for its next report. If the
     user chooses to stop, or nobody can answer (an unattended run), end the
     run and leave the branch as it is.
   - `LANDED`: verify it (next step).
   - `STOPPED`: end the run. Do not skip to a later ticket.
4. **Verify a landed ticket** yourself:
   - `git log -1 --format=%s main` starts with `<ID>: `, it is the only new
     commit since the previous ticket, and its message has a
     `World-breaking:` line;
   - the ticket file on `main` says `- Status: done` and its Result section is
     filled in;
   - the working tree is clean and no `ticket/<ID>-*` branch is left.

   If anything fails, end the run and report it; never repair a landed commit
   by rewriting history.
5. **Record.** Complete the ticket's task, append the outcome to the run log,
   and send the user a short progress line with SendUserMessage: ticket,
   commit, the key numbers from the summary, review rounds, follow-up tickets.
   New follow-up tickets that fall inside the target set join the plan (add
   tasks for them).
6. **Milestone.** If `docs/roadmap.md` now marks a milestone as met, send a
   short milestone summary (its tickets, the exit-criterion evidence, images
   worth opening) and continue.
7. **Target.** If the target is reached, the run ends.

Tickets run strictly one at a time: each builds on the previous commit on
`main`, and each must pass `mise run ci`, the same checks CI runs on Linux.

## Phase 3 — Finish

Report in Japanese: the tickets that landed in this run with their commits on
the local `main` (not pushed), why the run ended (target reached, waiting for
the owner, stopped, nothing left), the open question or open findings if any,
the next ticket, and the run log path. Leave tasks that were not reached
pending.

## Fallback without nested subagents

If `ticket-implementer` cannot be started, or it reports that it cannot start
the `ticket-reviewer` subagent (nesting may be disabled with
`CLAUDE_CODE_MAX_SUBAGENT_SPAWN_DEPTH=1`), tell the user and continue in this
conversation instead: for each ticket, invoke the `implement-ticket` skill with
the Skill tool and the ticket ID, and keep every other step of Phase 2. The
reviewer is then started directly from this conversation. The run log keeps
the run on track when the conversation is compacted.

## Never

- Implement, refine or review a ticket in this conversation, except through the
  fallback above.
- Answer an owner question yourself, approve a draft plan on the owner's
  behalf, or treat silence as approval.
- Skip a ticket that stopped and go on with later ones, or run two tickets at
  once.
- Relax any rule of `implement-ticket` to keep the run going: no landing
  without an approving review or a passing `mise run ci`, no golden value
  changed without a `World-breaking:` reason, no history rewriting.
- Run `git push`, `git pull` or `git fetch`, or add or change remotes.
- Edit tracked files yourself; the only file you write is the run log under
  `out/`.
