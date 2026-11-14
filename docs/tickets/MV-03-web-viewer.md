# MV-03 Viewer in the browser

- Status: done
- Milestone: MV
- Depends on: MV-02, MV-04, MV-06
- Crates: kuniumi-viewer, kuniumi-io
- Design: 実装ロードマップ / 技術スタック
- Decisions: D-10, D-11, D-19, D-23

## Goal

The same viewer runs in a web browser, so worlds can be inspected without
installing anything.

## Scope

In:
- Build for `wasm32-unknown-unknown` with WebGPU and a WebGL2 fallback through
  `wgpu`, plus the web platform glue it needs (`wasm-bindgen` family), listed
  in the D-10 table.
- Loading a world the user selects (file picker or drag and drop of the world
  directory); no server-side processing.
- A documented build script and a minimal HTML page; build output git-ignored.

Out:
- Hosting the page publicly (the repository is private until 1.0, D-03).
- Archives of a world directory (`.zip`, `.tar`): the directory itself is
  picked or dropped. Start options in the page's address (`--layer`,
  `--relief`, …) other than the GPU backend; screenshots and the benchmark
  in the browser.
- Opening a second world in the same page (winit allows one event loop per
  page): the page is reloaded instead.

## Implementation notes

Approved by the owner as refined: one ticket, and `web-serve` with the
pinned Python's `http.server`, recorded in D-11.

### What was checked while refining (scratch copy of the workspace)

- The library crates (`kuniumi-io`, `kuniumi-sim` and below) and today's
  viewer already compile for `wasm32-unknown-unknown` with the pinned 1.97.0
  toolchain. `std::fs` compiles there but every call fails, and
  `std::time::Instant::now` panics, so the web build reads no file and calls
  no `Instant`.
- `wasm-bindgen`, `wasm-bindgen-futures`, `js-sys` and `web-sys` are already
  in `Cargo.lock` (wgpu's and winit's dependencies on wasm32, 0.2.130,
  0.4.80, 0.3.107, 0.3.107). wgpu's `webgl` feature, needed for the WebGL2
  fallback, adds one package (`wgpu-core-deps-wasm`); `cargo deny check`
  passes (advisories, bans, licences, sources ok; duplicate-version warnings
  only, as today). Neither feature touches the native build.
- `wasm-bindgen-cli` 0.2.130 installs through mise's cargo backend in 21 s
  (`"cargo:wasm-bindgen-cli" = "0.2.130"`), and runs on the release `.wasm`
  of the binary (no `cdylib` needed; `main` becomes the start function).
  A release build of the viewer for wasm32 takes 23 s cold; after
  `wasm-bindgen` the module is 1.9 MB.
- `rust = { …, targets = ["wasm32-unknown-unknown"] }` in `mise.toml` makes
  `mise install` add the target, also to an installed toolchain that lacks
  it.
- On this machine, headless Chrome 157 has WebGPU only with
  `--enable-unsafe-webgpu --enable-features=Vulkan`; headless Firefox 158
  has none, so it takes the WebGL2 fallback. (Headed browsers are checked
  in the evidence.)

### Build tool and tasks (`mise.toml`, D-23)

`wasm-bindgen-cli`, not trunk: one tool that does one step, pinned in
`mise.toml` like cargo-deny; trunk would download its own `wasm-bindgen` at
build time and needs its own page conventions, and mise's registry name
`trunk` is another program. The CLI must equal the `wasm-bindgen` crate's
version, so the crate is pinned exactly (`=0.2.130`) in the viewer, and the
D-10 table says so (the only exact pin off the generation path, for this
reason).

- `[tools]`: `rust` gains `targets = ["wasm32-unknown-unknown"]`;
  `"cargo:wasm-bindgen-cli" = "0.2.130"`.
- `web`: `cargo build --release --locked -p kuniumi-viewer --target
  wasm32-unknown-unknown`, then `wasm-bindgen --target web --no-typescript
  --out-dir target/web target/wasm32-unknown-unknown/release/kuniumi-viewer.wasm`,
  then copy `crates/kuniumi-viewer/web/index.html` into `target/web/`. The
  output lies under `target/`, already git-ignored.
- `web-serve` (depends on `web`): serves `target/web` on
  `http://127.0.0.1:8000/` with the pinned Python's standard library,
  `python -m http.server --bind 127.0.0.1 --directory target/web 8000` (no
  script, no package; D-11 gains this development use). A browser does not
  load a WebAssembly module from a `file://` page, so some local server is
  needed; this one adds no tool (chosen by the owner).
- `web-check`, part of `check` (so of `ci`): `cargo clippy -p
  kuniumi-viewer --target wasm32-unknown-unknown --locked -- -D warnings`
  and `cargo doc -p kuniumi-viewer --no-deps --target
  wasm32-unknown-unknown --locked` with `RUSTDOCFLAGS=-D warnings` (the
  binary only: the tests run natively). The web code is compiled only for
  wasm32, so without this task nothing would lint it.

### Dependencies (`crates/kuniumi-viewer/Cargo.toml`)

Under `[target.'cfg(target_arch = "wasm32")'.dependencies]`, so the native
build is unchanged: `wgpu = { version = "30.0.2", features = ["webgl"] }`
(the same wgpu with the fallback added), `wasm-bindgen = "=0.2.130"`,
`wasm-bindgen-futures = "0.4.80"`, `js-sys = "0.3.107"`, `web-sys =
"0.3.107"` with the features `Document`, `HtmlCanvasElement`, `Performance`,
`Window` and `console`. No panic-hook crate (a three-line
`std::panic::set_hook` writes to the console) and no `web-time` (the page's
`performance.now()` through `web-sys`).

### Reading a world from memory (`kuniumi-io`)

```rust
/// Reads the world in the directory `dir` as [`read()`] does, opening each
/// file through `open` (given the file's path, `dir` joined with its name).
pub fn read_with<R: Read>(
    dir: &Path,
    open: impl FnMut(&Path) -> io::Result<R>,
) -> Result<World, ReadError>;
```

`read(dir)` becomes `read_with(dir, |path| File::open(path).map(BufReader::new))`;
`world.json` is read to its end from `open`, `read_binary` takes the reader.
The files asked for, their order and every error stay as today (an error of
`open` is `ReadError::Io` with the path). The browser hands each file's
bytes to `open` as a `Cursor<Vec<u8>>`.

### The page (`crates/kuniumi-viewer/web/index.html`, new)

A full-window `<canvas>`, a bar with a button for the directory picker
(`<input type="file" webkitdirectory>`), a status line and the keys and
mouse help; the whole page accepts a dropped directory (walked with
`webkitGetAsEntry`, or the dropped files themselves). JavaScript (about 100
lines, in the page) keeps only `world.json`, `grid.kgrid` and `*.kfield`
(the maps are not needed), reads them into `ArrayBuffer`s, and calls the
module's `open_world`. After a world is open, the picker is hidden;
another world needs a reload. `?backend=webgl` or `?backend=webgpu` in the
address forces one backend (default: WebGPU where the browser has an
adapter, else WebGL2).

### The web entry (`src/web.rs`, new, wasm32 only)

```rust
#[wasm_bindgen]
pub async fn open_world(
    canvas: HtmlCanvasElement,
    paths: Vec<String>,        // relative paths of the chosen files
    files: Vec<js_sys::Uint8Array>,
    backend: String,           // "", "webgl" or "webgpu"
    chosen_ms: f64,            // performance.now() when the files were chosen
) -> Result<String, JsValue>;   // the status line, or the error's message
```

1. `page::world_files(&paths)` finds the world; `world_dir::read_with`
   reads it, `open` copying one `Uint8Array` at a time into wasm memory, so
   the files' bytes are not all held twice.
2. `Scene::new` (mesh, layers, picker: the native code), `start_view` with
   the defaults of `args::parse([dir])`, `relief_note` to the console.
3. GPU, asynchronously (no `block_on` in a browser): the instance from
   `wgpu::util::new_instance_with_webgpu_detection` with the backends of
   `page::backends(backend)`; the surface from the canvas
   (`SurfaceTarget::Canvas`) before the adapter, which WebGL2 needs
   (`compatible_surface`); `gpu::request_device` (now `async`); the surface
   configured by `gpu::surface_config` (shared with the window, sizes
   clamped to `max_texture_dimension_2d`); `Renderer::new` unchanged.
4. winit's event loop with `EventLoopExtWebSys::spawn_app` (returns at
   once); the window on the same canvas
   (`WindowAttributesExtWebSys::with_canvas`, focusable, default actions
   prevented so keys and the wheel do not scroll the page).
5. Times (`performance.now()`, ms, from `chosen_ms`): files read by the
   browser, world read, scene built, GPU ready, first frame presented; and
   the wasm memory's size after loading (`wasm_bindgen::memory()`), printed
   to the console and returned in the status line with the backend
   (`BrowserWebGpu` or `Gl`) and the adapter's name.

A panic is written to the console (`std::panic::set_hook`) and shown in
the status line.

### Memory for 650,000 cells

Estimate (n = 650,000): the world 170 MB (the grid 101 MB, 13 fields
68 MB), the mesh 122 MB, the picker about 8 MB, one layer's vertex colours
18 MB, and while reading the largest file's bytes once more (101 MB); the
GPU upload's byte copies of the vertices and indices (120 MB) are passed
straight to the GPU. Peak wasm memory about 0.45–0.6 GB against wasm32's
4 GiB; the JavaScript side holds the chosen files (170 MB) until
`open_world` returns. Measured and recorded in the evidence.

### Shared code (`app.rs`, `gpu.rs`, `main.rs`)

- `app.rs`: `App` takes an optional `Prepared { surface, config,
  view_format, renderer }`: native `resumed` creates the window and the
  GPU as today (`block_on`); web `resumed` creates the window on the canvas
  and attaches the prepared GPU. `Action::Quit` (Q, Esc) does nothing in
  the browser. The title goes to `document.title` in the browser (winit
  sets only the canvas's `alt`). Lines for the terminal (legend, relief,
  the clicked cell) go through `say`, the console in the browser. The
  benchmark stays native.
- `gpu.rs`: `block_on`, `Offscreen` and the screenshots are native only;
  `request_device` is `async`; `surface_config` moves here from `app.rs`.
- `main.rs`: `run`, the screenshot, the benchmark and `crosshair` native
  only; `fn main() {}` for wasm32 (wasm-bindgen calls it at start); `say`.
- `page.rs` (new; compiled for wasm32 and for tests, so native tests cover
  it):
  - `world_files(paths: &[String]) -> Result<WorldFiles, String>`: the one
    `world.json` among the chosen paths (separator `/`, as browsers give
    `webkitRelativePath`; a bare name for dropped files), the directory
    that holds it (its name for the scene, `world` if the files were chosen
    without one) and its files by name; files in other directories ignored.
    Errors: no `world.json` ("choose a world directory written by
    `kuniumi generate`"), several ("… hold N worlds: …; choose one").
  - `backends(name: &str) -> Result<wgpu::Backends, String>`: `""` →
    `BROWSER_WEBGPU | GL`, `webgpu`, `webgl`, else an error.

### Tests (native, `mise run ci`)

- `kuniumi-io` (`read_tests.rs`): a written 1,000-cell world read through
  `read_with` from an in-memory map of its files equals `read()` bit for
  bit (every field's and the grid's bits, the metadata); `open` is asked
  for `world.json`, `grid.kgrid` and each listed `.kfield`, in that order,
  and for nothing else; a missing file gives `ReadError::Io` naming its
  path, a truncated `.kfield` `ReadError::Invalid`, as from files.
- `page`: `world_files` on a picked directory (`650k/world.json`, …,
  `650k/plate_id.png`), on a nested pick (`out/MV-03/650k/…` with a sibling
  directory without `world.json`), on bare names, with no and with two
  `world.json`; `backends` for each name and an unknown one.
- The existing tests keep passing unchanged (native code paths).

The browser code itself (`web.rs`, the page) cannot run in `cargo test`;
`web-check` lints it and the evidence exercises it in both browsers.

### Size

About 950 changed lines of code, tests and the page (viewer 600, page 170,
`kuniumi-io` 140, `mise.toml` and `Cargo.toml` 50) and 350 of documents.
Over the guideline of 800: a split (the reader of `kuniumi-io` as a ticket
of its own) was offered; the owner kept one ticket, as every part serves
the one feature.

### Documents

`crates/kuniumi-viewer/README.md` (new): native use in short, the browser
build (`mise install`, `mise run web`, `mise run web-serve`), browsers
(WebGPU in Chromium, the WebGL2 fallback, `?backend=`), choosing a world,
the controls and what differs from the window (Q and Esc, the console),
memory. `README.md` (the browser in one line), `docs/decisions.md` (D-10
row: `wasm-bindgen` pinned exactly to match the CLI; D-23: the wasm32
target, `wasm-bindgen-cli`, the web tasks; D-11: serving the page with
the pinned Python), `CHANGELOG.md` (Added), `docs/roadmap.md` if the MV
criterion is met.

## Acceptance criteria

- [x] Chromium-based browser (Chrome, WebGPU): `mise run web-serve`, then
      the directory `out/MV-03/650k` (`kuniumi generate --seed 1 --cells
      650000`) chosen with the picker shows the globe; the status line
      names the backend `BrowserWebGpu`; the load time from the choice to
      the first frame is recorded with its parts (files read, world read,
      scene, GPU, first frame) and the wasm memory's size after loading. A
      dropped directory loads as well.
- [x] Firefox (WebGL2 fallback): the same world loads with the backend `Gl`
      (by default where Firefox has no WebGPU, else with `?backend=webgl`);
      load time and memory recorded. Chrome with `?backend=webgl` recorded
      too.
- [x] The same layers and inspection as the native viewer: the page runs
      the native `Scene`, `Layers`, `press`, `Picker` and `panel` code (no
      second implementation). In both browsers, screenshots of the page:
      `surface_height` with relief and sea; another layer after Right; the
      panel for the cell under the pointer moved to a pixel, whose cell id
      and values equal what `kuniumi-viewer --screenshot --cursor` prints
      for the same view, size and pixel. Up/Down, S, P and R act as in the
      window, a click writes the cell's text to the console, Q and Esc do
      nothing.
- [x] `world_dir::read_with` and `page` tests pass (see Tests).
- [x] The build instructions in the crate README work from a clean
      checkout: in a fresh `git clone` of the branch, `mise install`,
      `mise run web` and `mise run web-serve` serve the page, which loads
      the world.
- [x] `cargo tree -e normal -p kuniumi-sim -p kuniumi-io -p kuniumi-cli
      --target all` shows no `wgpu`, `winit`, `wasm-bindgen`, `js-sys` or
      `web-sys`.
- [x] `mise run ci` passes, with `web-check` part of it.

## Verification

```sh
mise install
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 650000 --out out/MV-03/650k
mise run web-serve      # then open http://127.0.0.1:8000/ and choose out/MV-03/650k
# http://127.0.0.1:8000/?backend=webgl forces the WebGL2 fallback
cargo run --release -p kuniumi-viewer -- out/MV-03/650k --cursor 640,400 --screenshot out/MV-03/native.png
```

Screenshots of the page are taken through WebDriver (Firefox, geckodriver)
and the DevTools protocol (Chrome), with the window's size and the pointer
set, and compared with the native screenshot of the same view.

## Result

What was done:
- `kuniumi-io`: `world_dir::read_with(dir, open)`; `read(dir)` is
  `read_with` on the file system, `read_binary` takes the opened reader.
  Tests `reads_through_any_opener`, `opener_errors_name_their_file`,
  `spoiled_files_read_from_memory_fail_as_from_disk`.
- `kuniumi-viewer`: `page.rs` (new; `world_files`, `WorldFiles::{index,
  name, wanted}`, `backends`, `LoadTimes::status`, with 5 tests), `web.rs`
  (new; `wanted_files`, `open_world`), `app.rs` (`Prepared`, `prepare`,
  `spawn`, `App::new`, `open` with the prepared GPU or the native one,
  `fail`, Q/Esc ignored in a browser, `document.title`, `say`, the first
  frame's callback, resizes clamped to the largest texture side), `gpu.rs`
  (`request_device` async, `surface_config` moved here, `block_on` and
  `Offscreen` native only), `main.rs` (`say`, `fn main() {}` for wasm32,
  native-only `run` and `crosshair`), `args.rs` (`USAGE` native only;
  `default_for`, the page's default view, review R1-1),
  `web/index.html` (new), `README.md` (new), `Cargo.toml` (wasm32
  dependencies; `Cargo.lock` gains `wgpu-core-deps-wasm` only).
- `mise.toml`: the wasm32 target, `cargo:wasm-bindgen-cli` 0.2.130, the
  tasks `web`, `web-serve`, `web-check` (in `check`).
- Documents: D-10 (row), D-11 (serving), D-23 (target, tool, tasks),
  `design-ja.md` (技術スタック: 開発環境, ビューア), `README.md`, `CLAUDE.md`
  (workspace block: the browser build is no longer planned), `CHANGELOG.md`
  (Added), `docs/conventions.md` (the one use of Python outside `tools/`),
  `docs/roadmap.md` (MV met: with this ticket a 650,000-cell world can be
  explored natively and in the browser; MV-07, a draft outside the row's
  MV-01 … MV-06, refines the panel).

Deviations from the plan:
- `wanted_files` (exported) lets the page read only the files the reader
  may ask for, decided in Rust (`WorldFiles::wanted`) rather than by a
  name filter in JavaScript; `open_world` takes a `report` function for
  the status line after the first frame and for errors that stop the
  window; the load times and their line are `page::LoadTimes` (tested)
  rather than formatted in `web.rs`.
- `web-sys` also needs the features `Element` and `HtmlElement` (the
  canvas's CSS size and `focus`) and `Node` (the status line written by the
  panic hook, R1-2).
- The adapter's name is left out of the status line when the browser
  hides it (Chrome's WebGPU gives none).
- Added after review round 1: the panic hook also writes the panic into
  the page's status line, and the page shows any other error that stops
  the module (`error` and `unhandledrejection` listeners), as a panic
  aborts the module and its promise never settles (R1-2); keys with Ctrl,
  Alt or Meta and the function keys are kept from the canvas, whose winit
  listener prevents every key's default action, so reload, find and zoom
  stay the browser's (R1-5).
- Added after review round 2: the panic hook is installed by the wasm32
  `main`, when the module starts, so it covers `wanted_files` too (R2-2);
  the page keeps a panic's message when the call that trapped throws, and
  replaces any other status, a recoverable error included, with "the
  viewer stopped"; an error while walking a dropped directory is shown as
  a load error, from which another drop recovers (R2-3).
- Size (before review): about 1,290 changed lines of code, tests, the page and the crate
  README (the tests 480 of them), and 300 of documents, against the
  estimate of 950 and 350.

Evidence (this machine: NVIDIA GeForce RTX 3060 Ti, Ubuntu GNOME on
Wayland; world `out/MV-03/650k` from `kuniumi generate --seed 1 --cells
650000`, 650,000 cells, 13 fields; a 20,000-cell world `out/MV-03/20k`
too; the pages driven through the DevTools protocol (Chrome) and WebDriver
(geckodriver, Firefox) by throwaway scripts, headed, at a 1280 × 800
viewport unless stated):

1. Chrome 157 (WebGPU, no flags): the directory chosen with the picker
   (DevTools `DOM.setFileInputFiles` on the `webkitdirectory` input):
   status `650k: 650000 cells, 13 fields; BrowserWebGpu; files 0.20 s,
   world 0.27 s, scene 0.20 s, GPU 0.63 s, first frame 1.32 s after the
   choice; memory 542 MB`. Over seven loads: first frame 1.03–1.47 s
   (files 0.19–0.21 s, world 0.27–0.28 s, scene 0.20 s, GPU 0.32–0.78 s),
   wasm memory 542 MB every time. The directory dropped on the page
   (DevTools `Input.dispatchDragEvent` with the directory's path, walked
   by `webkitGetAsEntry`): the same world, first frame 1.31 s
   (`out/MV-03/chrome_650k_dropped.png`).
2. Firefox 158 (no `navigator.gpu`, so the WebGL2 fallback by default):
   status `world: 650000 cells, 13 fields; Gl on NVIDIA GeForce GTX 980,
   or similar; files 0.17 s, world 0.27 s, scene 0.22 s, GPU 0.57 s,
   first frame 1.44 s after the choice; memory 589 MB` (a second load:
   1.24 s, 589 MB). geckodriver refuses a directory for the
   `webkitdirectory` input ("File not found"), so for Firefox the input
   was switched to `multiple` in the page and the directory's 30 files
   given; the page then sees bare names (world `world`), the case of
   dropped files. Chrome with `?backend=webgl`: `Gl on ANGLE (NVIDIA
   Corporation, NVIDIA GeForce RTX 3060 Ti/PCIe/SSE2, OpenGL ES 3.2)`,
   first frame 1.31 s, memory 589 MB.
3. Same layers and inspection: the page runs `Scene`, `Layers`, `press`,
   `Picker`, `panel` and the window code of `app.rs`; only the GPU's
   preparation and the event loop's start differ (`web.rs`).
   - The pointer moved to (640.5, 400.5) and clicked: the console holds
     `cell 214129 of 650000, lat 19.9256 deg, lon 0.0805 deg`, the 13
     field lines and the layer line, identical line for line to what
     `kuniumi-viewer out/MV-03/650k --cursor 640,400 --screenshot
     out/MV-03/native_650k_cursor.png` prints; in Firefox (WebGL2) the
     same lines, and in Chrome with `?backend=webgl` the same cell.
   - Screenshots against the native screenshot of the same view and
     pixel (`out/MV-03/native_650k_cursor.png`; compared above the page's
     status bar, the crosshair's box left out): Firefox
     (`firefox_650k_cursor.png`) 95.5 % of the pixels identical, all but a
     few within 2 levels (largest difference 19, at multisampled edges);
     after Right × 3, `plate_id` (`firefox_650k_plate_id.png` against
     `native_650k_plate_id_cursor.png`): 96.0 % identical, largest 56.
     Chrome (`chrome_srgb_cursor.png`, `chrome_srgb_plate_id.png`): every
     pixel within 8 levels and none identical, because Chrome's capture
     transforms colours: a plain CSS `#303030` square is captured as 51,
     `rgb(34, 87, 159)` as (39, 88, 158), the same shifts as on the
     canvas (also with `--force-color-profile=srgb`).
   - Looked at: `chrome_650k_cursor.png` (the globe as natively: the
     collision belt in yellow-green across the lower left, the three
     continents, light-blue shelves, the panel with cell 214129 and the
     colour bar), `chrome_srgb_plate_id.png` and `firefox_650k_plate_id.png`
     (the plates in the categorical palette, tinted blue by the sea, the
     panel with `layer 3/13: plate_id`), `firefox_650k_nopanel_nosea.png`
     (P and S: no panel, no sea, the plates with their raised continents
     at relief ×50), `chrome_srgb_zoomed.png` and
     `firefox_650k_dragged_zoomed.png` (after the same drag and wheel in
     both browsers: the same view, cell 233327 in the panel),
     `chrome_20k.png` (hexagonal cells visible), `chrome_dpr2_real.png`
     (Chrome with `--force-device-scale-factor=2`: canvas 1854 × 874 for
     927 × 437 CSS px, sharp, the panel at text scale 2).
   - Keys (titles read from `document.title`): Up `relief x50`, Down ×2
     `relief x10`, S `sea off`, P hides the panel, R after a drag and a
     wheel restores the view (the screenshot after R equals the one before
     the drag, bit for bit), Q and Esc change nothing (the screenshot after
     them equals the one before, bit for bit; frames go on), Right/Left
     cycle `1/13 cell_area` … `3/13 plate_id` and back; the console lines
     of the legends and relief steps as in the terminal.
4. Tests: `kuniumi-io` 120 (3 new), `kuniumi-viewer` 68 (5 new in
   `page`). Test-first: the three `read_with` tests were run against an
   `unimplemented!()` stub and failed (the stub was replaced before the
   commit dd38033); the four `world_files`/`backends` tests failed against
   the stubs of commit 3228c96. `LoadTimes::status` was written with its
   test.
5. Clean checkout: `git clone --branch ticket/MV-03-web-viewer` into a
   scratch directory, `mise trust`, `mise install`, `mise run web` (23.7 s;
   `target/web/` with `index.html`, `kuniumi-viewer.js` 148 kB,
   `kuniumi-viewer_bg.wasm` 5.1 MB), `mise run web-serve`: the page served
   from the clone (its `.wasm` byte for byte) loaded the 650,000-cell world
   in Chrome, first frame 1.47 s (`out/MV-03/clean_clone_chrome.png`).
   `mise install` re-adds the wasm32 target when it is missing from the
   installed toolchain (checked by removing it while refining).
6. `cargo tree -e normal -p <crate> --target all`: no line with `wgpu`,
   `winit`, `wasm-bindgen`, `js-sys` or `web-sys` for `kuniumi-sim`,
   `kuniumi-io`, `kuniumi-cli`, `kuniumi-math`, `kuniumi-sphere`,
   `kuniumi-world`, `kuniumi-config`.
7. `mise run ci` passes (log `out/MV-03/ci.log`, after the round-1 fixes):
   fmt, clippy, doc, `web-check` (clippy and rustdoc of the viewer for
   wasm32, warnings as errors), 1,225 tests passed in debug and release
   together (10 ignored: the existing timing tests), MSRV 1.85.0 without
   the viewer, cargo-deny (advisories, bans, licences, sources ok;
   duplicate-version warnings only).
8. Round-1 fixes, checked in the browsers on the rebuilt page (Chrome
   with the cache disabled): a copy of the 20,000-cell world in a
   directory named `-old` loads (title `kuniumi-viewer - -old - 13/13
   surface_height …`; test `default_view_for_any_directory_name` for
   `-old`, `--help`, `--relief`, `-h` and an empty name, failing against a
   stub first). A throwaway build that panics in `open_world` on
   `?backend=panic` (not committed) shows `error: the viewer stopped:
   panicked at crates/kuniumi-viewer/src/web.rs:111:9: a test panic;
   reload the page` in the status line; an error thrown by a page script
   shows `error: the viewer stopped: Uncaught Error: a test error; …`.
   Keydown events dispatched on the canvas: before the fix F5 and Ctrl+R
   were default-prevented and Ctrl+Right changed the layer; after it, in
   Chrome and Firefox, F5 and Ctrl+R are not prevented, Ctrl+Right leaves
   the layer, Right is prevented and changes it. The Chrome interaction
   run repeated on the rebuilt page gives the same titles, cell 214129,
   and the screenshots after R and after Q/Esc equal the one before, bit
   for bit; Firefox: first frame 1.40 s, 589 MB, cell 214129.
9. Round-2 fixes, in Chrome: a throwaway build that panics in
   `wanted_files` when a chosen file is named `panic` (not committed)
   shows `error: the viewer stopped: panicked at
   crates/kuniumi-viewer/src/web.rs:75:9: a test panic in wanted_files;
   reload the page` (before the page's fix the synchronous trap replaced it
   with `error: RuntimeError: unreachable`); the real build, given the same
   file, says `error: no world.json among the 1 files chosen: …`. A drop
   whose directory walk fails (a synthetic drop event whose entry's
   `file` fails) shows `error: the dropped files cannot be read: Error:
   unreadable`; after that recoverable error and a picked directory
   without `world.json`, an error thrown by a page script shows `error:
   the viewer stopped: Uncaught Error: a test error; …`. The 650,000-cell
   world then loads as before (first frame 1.20–1.36 s, 542 MB, cell
   214129 on a click).

Notes:
- During the headed runs the desktop's own mouse pointer, resting over the
  browser windows, sends pointer events of its own (seen in a trace at
  device pixel ratio 2); one Firefox screenshot taken that way showed "no
  cell under the cursor" and was taken again.

Review:
- Round 1 (CHANGES_REQUESTED): the reviewer re-ran `mise run ci`, both
  browsers (Chrome WebGPU 1.12–1.58 s and 542 MB; Firefox 1.42 s and
  589 MB; Chrome WebGL2 1.37 s), the 15 cell lines against the native
  output, the screenshot comparisons (the same figures), the keys, the
  error paths and the clean clone, and the native benchmark (59.2 fps).
  R1-1 major, a world directory whose name starts with `-` failed in the
  page ("the default view cannot be made": the name went through the
  argument parser) (fixed: `args::default_for`, test). R1-2 minor, a panic
  wrote only to the console and the status line kept "Building the
  globe…" (fixed: the hook writes the status line; the page's `error`
  and `unhandledrejection` listeners). R1-3 minor, the roadmap's MV row
  was not set to met (fixed). R1-4 nit, `docs/conventions.md` and
  `CLAUDE.md` rule 3 still said Python only under `tools/` (fixed in
  `conventions.md`; `CLAUDE.md` left as it is: its rule cites D-11, which
  now records the exception, and `CLAUDE.md` is not edited on a review's
  request). R1-5 nit, winit's default prevention could swallow F5 and
  Ctrl+R, which the help tells users to use (confirmed: they were
  prevented; fixed).
- Round 2 (APPROVE): the reviewer re-ran `mise run ci` and both browsers
  (Chrome WebGPU 1.12–1.57 s, 542 MB; Firefox 1.34 s, 589 MB; Chrome
  WebGL2 1.16 s), with its own panicking builds (in `open_world` and in
  the event loop) and directories named `-old` and `--help`; R1-1 to R1-5
  resolved (R1-4's rebuttal for `CLAUDE.md` accepted). Nits, all fixed:
  R2-1, the deviation about `web-sys` features left out `Node`; R2-2, the
  panic hook was installed only by `open_world`, not for `wanted_files`;
  R2-3, the page kept a recoverable error over a later fatal one, and a
  failed directory walk of a drop was reported as the viewer stopping.
- Round 3 (APPROVE, no findings): the reviewer re-ran `mise run ci` on
  de4d70a, rebuilt the page from a clean viewer build (the same `.wasm`
  byte for byte), both browsers (Chrome WebGPU 1.37 s and 542 MB, drop
  1.40 s; Firefox 1.39 s and 589 MB; Chrome WebGL2 1.34 s), the 15 cell
  lines against the native output, the screenshot comparisons, its own
  panicking builds (in `wanted_files` and in `open_world`) and the page's
  error paths; R2-1 to R2-3 resolved.

Follow-ups: none.
