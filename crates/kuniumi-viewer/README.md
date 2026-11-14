# kuniumi-viewer

A world directory written by `kuniumi generate` as an interactive 3D globe,
in a window or in a web browser. Every field of the world is a layer, in the
colours of its PNG map; the relief is exaggerated over a translucent sea
surface; the cell under the pointer is highlighted and its values are shown
in a panel with the layer's legend.

## In a window

Rust 1.87 or later and a GPU with Vulkan, Metal or Direct3D 12:

```sh
cargo run --release -p kuniumi-cli -- generate --seed 1 --cells 650000 --out out/world
cargo run --release -p kuniumi-viewer -- out/world
cargo run --release -p kuniumi-viewer -- --help    # start view, screenshots, benchmark
```

## In a browser

The same viewer, built for WebAssembly, runs in a page that opens a world
directory chosen on the user's machine. Nothing is uploaded: the browser
reads the files and the page draws them.

### Build and serve

From a clone of the repository, with [mise](https://mise.jdx.dev):

```sh
mise trust            # once in a fresh clone: let mise use its mise.toml
mise install          # the pinned Rust with its wasm32 target, wasm-bindgen-cli, Python
mise run web          # builds the page into target/web/
mise run web-serve    # builds it and serves target/web/ on http://127.0.0.1:8000/
```

`mise run web` compiles the viewer for `wasm32-unknown-unknown` in release
mode, generates the JavaScript glue with `wasm-bindgen` (its version is
pinned in `mise.toml` and must equal the `wasm-bindgen` crate's) and copies
`web/index.html` next to it. `target/web/` holds three files
(`index.html`, `kuniumi-viewer.js`, `kuniumi-viewer_bg.wasm`, about 5 MB)
and can be served by any static file server; `mise run web-serve` uses the
standard library of the pinned Python and serves only on this machine. A
page opened as a `file://` address does not work: browsers do not load a
WebAssembly module from it.

### Open a world

Open <http://127.0.0.1:8000/>, then either press **Choose a world
directory** and pick the directory (for example `out/world`, the one that
holds `world.json`), or drop the directory onto the page. The page reads
`world.json`, `grid.kgrid` and the `.kfield` files (not the PNG maps), and
the status line at the bottom then shows the GPU backend and the load
times. One world per page: reload the page to open another.

### Browsers

- WebGPU where the browser offers it (Chromium-based browsers; tested with
  Chrome 157 on Linux). Elsewhere the page falls back to WebGL2 (tested with
  Firefox 158 on Linux, which has no WebGPU there).
- `?backend=webgl` or `?backend=webgpu` at the end of the address forces one
  of them (`http://127.0.0.1:8000/?backend=webgl`).
- A world of 650,000 cells takes about 0.55 GB of the page's WebAssembly
  memory (WebGPU; 0.59 GB with WebGL2) of the 4 GiB a 32-bit module can
  have, and loads in about 1.0–1.5 s on a desktop machine (ticket MV-03).

### Controls

As in the window: drag to turn the globe and use the wheel to zoom;
Left/Right change the layer, Up/Down the relief, S shows or hides the sea,
P the panel, R resets the view. Differences from the window: Esc and Q do
nothing (close the tab instead); a click writes the cell's values to the
browser's console (the developer tools) instead of the terminal, as are the
legends and the relief steps; the title of the window is the title of the
page.
