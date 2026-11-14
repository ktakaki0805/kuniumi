# MV-07 Panel height bounded by the window

- Status: draft
- Milestone: MV
- Depends on: MV-06
- Crates: kuniumi-viewer
- Design: 実装ロードマップ / 技術スタック
- Decisions: D-19

## Goal

The panel of `kuniumi-viewer` (MV-06) fits the window's height however many
fields a world has, so its legend and colour bar stay visible.

## Scope

In:
- A bound on the panel's height from the viewport's height and the text
  scale: when the lines do not fit, the field lines are cut (with a line
  saying how many are hidden) or shown in columns, the header, the layer
  line and the colour bar always kept.
- The choice between cutting and columns, and whether the field shown
  first (the layer's) is always among the visible lines.

Out:
- Scrolling or selecting fields with the mouse; a panel that moves away
  from the cursor.

## Implementation notes

Found in MV-06's review (R1-5): only the panel's width follows the viewport.
With the 13 fields of the M2 worlds it is 208 px tall at scale 1 and 416 px
at scale 2; later milestones add many fields (climate, hydrology, biomes),
and in a short window the bar and its labels are cut off at the bottom. `P`
hides the panel meanwhile. To be refined: the rule for which lines are kept,
tested without a GPU as `panel::layout` is.

## Acceptance criteria

- [ ] For every viewport height down to a stated minimum, the panel's
      quads lie inside the viewport and the layer line, the colour bar and
      its labels are among them (tests on a world with many fields).
- [ ] A screenshot of a world with many fields in a short window is saved
      under `out/MV-07/` and described.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-viewer -- out/MV-02/650k --size 1280x200 --cursor 640,100 --screenshot out/MV-07/short.png
```

## Result

_(filled in when the ticket is implemented)_
