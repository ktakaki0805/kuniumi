# MV-03 Viewer in the browser

- Status: draft
- Milestone: MV
- Depends on: MV-02, MV-04, MV-06
- Crates: kuniumi-viewer
- Design: 実装ロードマップ / 技術スタック
- Decisions: D-10, D-19

## Goal

The same viewer runs in a web browser, so worlds can be inspected without
installing anything.

## Scope

In:
- Build for `wasm32-unknown-unknown` with WebGPU and a WebGL2 fallback through
  `wgpu`, plus the web platform glue it needs (`wasm-bindgen` family), listed
  in the D-10 table.
- Loading a world the user selects (file picker or drag and drop of the world
  directory or an archive of it); no server-side processing.
- A documented build script and a minimal HTML page; build output git-ignored.

Out:
- Hosting the page publicly (the repository is private until 1.0, D-03).

## Implementation notes

To be refined: build tool (`wasm-bindgen-cli` or `trunk`), pinned in
`mise.toml` together with the `wasm32-unknown-unknown` target of the Rust tool
and a `mise run` task for the web build (D-23); how a directory is handed to
the browser; memory limits for 650,000-cell worlds.

## Acceptance criteria

- [ ] The page loads a 650,000-cell world in a current Chromium-based browser
      (WebGPU) and in Firefox (WebGL2 fallback); load time recorded.
- [ ] The same layers and inspection as the native viewer.
- [ ] The build instructions in the crate README work from a clean checkout.
- [ ] `mise run ci` passes.

## Verification

Follow the build instructions in the crate README and open the page in both
browsers.

## Result

_(filled in when the ticket is implemented)_
