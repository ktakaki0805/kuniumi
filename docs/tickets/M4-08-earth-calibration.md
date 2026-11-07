# M4-08 Earth calibration harness

- Status: draft
- Milestone: M4
- Depends on: M4-07
- Crates: kuniumi-cli, `tools/` (Python)
- Design: 気候モデル / 較正：地球を正解として使う; 検証方法
- Decisions: D-11, D-16

## Goal

Run the climate model on Earth's real topography and score it against observed
climate, so parameters are tuned against reality rather than by eye.

## Scope

In:
- Datasets: ETOPO 2022 topography and the Köppen–Geiger maps of Beck et al.
  (2018); choose observed temperature and precipitation datasets here. Confirm
  and record each dataset's licence and citation in D-16 before use. Data is
  downloaded by a script into a git-ignored folder and never committed.
- Conversion: Python scripts under `tools/` (packages pinned in
  `tools/requirements.txt`, D-11) read the source formats (NetCDF, GeoTIFF)
  and write simple little-endian arrays with a small header for the Rust side.
- Import: resample Earth topography onto the grid.
- Scores: zonal-mean temperature and precipitation errors (the Köppen agreement
  score is added by M5-04).
- A short report (text and images) produced by one command; parameters adjusted
  and the improvement recorded; `docs/roadmap.md`: M4 marked as met.

## Implementation notes

To be refined: the observed climate datasets (ask the owner before adding any
download), the resampling method, score definitions.

## Acceptance criteria

- [ ] One command produces the calibration report from the downloaded data.
- [ ] Baseline scores and the scores after tuning recorded in the Result; zonal
      means broadly match observations (the M4 exit criterion).
- [ ] No data files committed; licences recorded in D-16.
- [ ] `mise run ci` passes.

## Verification

Run the documented download, conversion and report commands; view the report
images.

## Result

_(filled in when the ticket is implemented)_
