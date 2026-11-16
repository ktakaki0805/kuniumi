# M4-07 Climate stage and monthly outputs

- Status: draft
- Milestone: M4
- Depends on: M4-06, M4-09
- Crates: kuniumi-sim, kuniumi-io, kuniumi-cli
- Design: 気候モデル / 計算の流れ
- Decisions: D-22

## Goal

The Climate stage runs in the pipeline before Erosion, and erosion uses the
modelled precipitation instead of the placeholder.

## Scope

In:
- Stage wiring with the monthly loop (insolation → temperature → pressure and
  winds → currents (annual) → sea-surface temperature → moisture and
  precipitation) until it is periodic.
- Fields: `temperature_m01` … `temperature_m12` (K), `precipitation_m01` …
  `precipitation_m12` (m/s), wind, currents, sea-surface temperature, plus
  annual means for quick viewing.
- Erosion discharge from annual-mean precipitation.
- Rendering of the annual means and of two contrasting months.
- A planet whose energy balance has not converged within its 100 years
  (`energy_balance::Temperature::converged`, M4-02) keeps the last year's
  temperatures and is reported as a warning, with the last change; the
  warning says that the climate may be in transition (near an ice–albedo
  bifurcation the last year can lie mid-way into a snowball: Earth at
  1,296 W/m² changed by 11.5 K in its 100th year). Decided with the owner
  in M4-02. The temperatures have no fixed bound at the corners of the v1
  range (monthly means from about 95 K to about 11,600 K, M4-02).

## Implementation notes

To be refined: the periodicity criterion and maximum number of model years,
which fields are kept monthly versus annual only (output size), whether wind
and currents are written as two scalar fields each. The moisture of M4-06
(`moisture::monthly_moisture`) also gives monthly evaporation and column
vapour, and reports `converged`, `unsettled_steps` and `failed_solves`, to be
reported like the energy balance's convergence.

## Acceptance criteria

- [ ] `kuniumi generate` writes all climate fields; images viewed and described.
- [ ] Golden fingerprints for 10,000 cells, seed 1 (intentional world-breaking
      update of the erosion outputs, recorded in the changelog).
- [ ] Timing for 650,000 cells recorded.
- [ ] `mise run ci` passes.

## Verification

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --out out/M4-07
```

## Result

_(filled in when the ticket is implemented)_
