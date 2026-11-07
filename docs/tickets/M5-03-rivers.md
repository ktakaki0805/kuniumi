# M5-03 Rivers, seasonal rivers and vector export

- Status: draft
- Milestone: M5
- Depends on: M5-02
- Crates: kuniumi-sim (stages/hydrology), kuniumi-io
- Design: 水文（川・湖）; 高解像度化と出力 / ゲームへの出力
- Decisions: none

## Goal

Rivers as data a game can use: which cells carry a river, whether it flows all
year, and polylines with discharge for drawing.

## Scope

In:
- River cells: annual discharge above a threshold (documented, scaled by cell
  size).
- Seasonal rivers (wadis): river cells whose monthly runoff is positive only in
  some months.
- Polylines from sources to outlets or lakes, split at confluences, with
  discharge per vertex; drawn width ∝ √discharge.
- Export: `rivers.geojson` with longitude and latitude in degrees (the export
  boundary; internal values stay in radians), properties for discharge and
  seasonality; `rivers.png` updated to use the polylines.

Out:
- Projections other than longitude–latitude (M8-03).

## Implementation notes

To be refined: threshold, polyline simplification, handling of rivers crossing
the antimeridian in GeoJSON.

## Acceptance criteria

- [ ] Every polyline ends at the ocean, a lake, or a confluence; none crosses
      another.
- [ ] `rivers.geojson` is valid GeoJSON (checked by a test that parses it back
      with `serde_json`).
- [ ] River map viewed and described (seasonal rivers in dry regions).
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M5-03
```

## Result

_(filled in when the ticket is implemented)_
