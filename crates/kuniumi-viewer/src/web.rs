//! The viewer in a browser (ticket MV-03): the functions the page
//! (`web/index.html`) calls. The page asks [`wanted_files`] which of the
//! files the user chose to read, reads them, and hands them to
//! [`open_world`], which reads the world from them
//! (`kuniumi_io::world_dir::read_with`), builds the same [`Scene`] as the
//! native viewer, prepares the GPU on the page's canvas (WebGPU, else
//! WebGL2) and starts the same window code ([`app::spawn`]) on it.
//!
//! Nothing here reads a file or the clock of `std` (neither works in a
//! browser): the files come from the page, and times from the page's
//! `performance.now()`.

use std::io::{self, Cursor};
use std::path::Path;

use wasm_bindgen::JsCast;
use wasm_bindgen::prelude::*;
use web_sys::HtmlCanvasElement;
use winit::platform::web::WindowAttributesExtWebSys;
use winit::window::Window;

use crate::args;
use crate::page::{self, LoadTimes};
use crate::{Scene, app, relief_note, say, start_view};

/// Milliseconds since the page started (`performance.now()`).
fn now_ms() -> f64 {
    web_sys::window()
        .and_then(|w| w.performance())
        .map_or(0.0, |p| p.now())
}

/// Sets the page's title (the window title of the native viewer).
pub fn set_title(title: &str) {
    if let Some(document) = web_sys::window().and_then(|w| w.document()) {
        document.set_title(title);
    }
}

/// Sends panics to the console and to the page's status line (the element
/// `status` of `web/index.html`). A panic aborts the module for good (no
/// unwinding on wasm32) and the page's promise never settles, so the page
/// learns of it only here. Called by `main` when the module starts, so it
/// covers every export.
pub fn install_panic_hook() {
    std::panic::set_hook(Box::new(|info| {
        web_sys::console::error_1(&format!("kuniumi-viewer panicked: {info}").into());
        let status = web_sys::window()
            .and_then(|w| w.document())
            .and_then(|d| d.get_element_by_id("status"));
        if let Some(status) = status {
            status.set_text_content(Some(&format!(
                "error: the viewer stopped: {info}; reload the page"
            )));
        }
    }));
}

/// Calls the page's function `report` with `text`.
fn tell(report: &js_sys::Function, text: &str) {
    // An exception of the page's own function is the page's to show.
    let _ = report.call1(&JsValue::NULL, &JsValue::from_str(text));
}

/// The indices among `paths` (the relative paths of the files the user
/// chose) of the files [`open_world`] needs ([`page::WorldFiles::wanted`]).
///
/// # Errors
///
/// The message of [`page::world_files`], if the files hold no world or
/// several.
#[wasm_bindgen]
pub fn wanted_files(paths: Vec<String>) -> Result<Vec<u32>, JsValue> {
    let world = page::world_files(&paths).map_err(|e| JsValue::from_str(&e))?;
    Ok(world.wanted().into_iter().map(|k| k as u32).collect())
}

/// Shows the world held by `files` (the bytes of the files at `paths`,
/// relative paths as the browser gives them) on `canvas`, with the GPU
/// backend `backend` of the page's address (`""`, `webgpu` or `webgl`,
/// [`page::backends`]). `chosen_ms` is the page's time when the user chose
/// the files; `report` is called with a line for the page's status: when
/// the first frame is presented (with the load times), and if the viewer
/// stops on an error. Returns the status line once the GPU is ready.
///
/// # Errors
///
/// The message of what failed: no world among the files, a file that is
/// not valid ([`kuniumi_io::world_dir::ReadError`]), no GPU.
#[wasm_bindgen]
pub async fn open_world(
    canvas: HtmlCanvasElement,
    paths: Vec<String>,
    files: Vec<js_sys::Uint8Array>,
    backend: String,
    chosen_ms: f64,
    report: js_sys::Function,
) -> Result<String, JsValue> {
    open(canvas, &paths, &files, &backend, chosen_ms, report)
        .await
        .map_err(|e| {
            say(&format!("kuniumi-viewer: {e}"));
            JsValue::from_str(&e)
        })
}

async fn open(
    canvas: HtmlCanvasElement,
    paths: &[String],
    files: &[js_sys::Uint8Array],
    backend: &str,
    chosen_ms: f64,
    report: js_sys::Function,
) -> Result<String, String> {
    let backends = page::backends(backend)?;
    let mut times = LoadTimes {
        files_ms: now_ms() - chosen_ms,
        ..LoadTimes::default()
    };
    let chosen = page::world_files(paths)?;
    if files.len() != paths.len() {
        return Err(format!(
            "{} files for {} paths from the page",
            files.len(),
            paths.len()
        ));
    }
    let dir = Path::new(chosen.name());
    let mark = now_ms();
    let world = kuniumi_io::world_dir::read_with(dir, |path| {
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        match chosen.index(name) {
            // One file at a time in the module's memory.
            Some(k) => Ok(Cursor::new(files[k].to_vec())),
            None => Err(io::Error::new(
                io::ErrorKind::NotFound,
                "not among the files chosen",
            )),
        }
    })
    .map_err(|e| e.to_string())?;
    times.read_ms = now_ms() - mark;

    let mark = now_ms();
    let args = args::default_for(dir.to_path_buf());
    let scene = Scene::new(dir.to_path_buf(), world, None)?;
    times.scene_ms = now_ms() - mark;
    say(&format!(
        "{}: {} cells, {} fields, mesh of {} triangles",
        chosen.name(),
        scene.world.grid.len(),
        scene.world.fields.len(),
        scene.mesh.indices.len() / 3
    ));
    let view = start_view(&args, &scene.relief);
    if let Some(note) = relief_note(args.relief, view.relief, &scene.relief) {
        say(&note);
    }

    let mark = now_ms();
    let instance = wgpu::util::new_instance_with_webgpu_detection(wgpu::InstanceDescriptor {
        backends,
        ..wgpu::InstanceDescriptor::new_without_display_handle()
    })
    .await;
    // Before the adapter, which WebGL2 finds through the surface.
    let surface = instance
        .create_surface(wgpu::SurfaceTarget::Canvas(canvas.clone()))
        .map_err(|e| format!("cannot draw into the page's canvas: {e}"))?;
    let adapter = instance
        .request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: Some(&surface),
            ..Default::default()
        })
        .await
        .map_err(|e| format!("no GPU adapter in this browser: {e}"))?;
    let info = adapter.get_info();
    let ratio = web_sys::window().map_or(1.0, |w| w.device_pixel_ratio());
    // The canvas's size in device pixels; winit reports changes later.
    let side = |css_px: i32| (f64::from(css_px.max(1)) * ratio).round() as u32;
    let (width_px, height_px) = (side(canvas.client_width()), side(canvas.client_height()));
    let prepared = app::prepare(&adapter, surface, width_px, height_px, &scene).await?;
    times.gpu_ms = now_ms() - mark;
    times.memory_bytes = wasm_bindgen::memory()
        .dyn_into::<js_sys::WebAssembly::Memory>()
        .ok()
        .and_then(|m| m.buffer().dyn_into::<js_sys::ArrayBuffer>().ok())
        .map_or(0, |b| u64::from(b.byte_length()));
    // Browsers may keep the adapter's name to themselves (WebGPU).
    let gpu = match info.name.as_str() {
        "" => format!("{:?}", info.backend),
        name => format!("{:?} on {name}", info.backend),
    };
    let (name, cells, fields) = (
        chosen.name().to_owned(),
        scene.world.grid.len(),
        scene.world.fields.len(),
    );
    let status = times.status(&name, cells, fields, &gpu);
    say(&status);

    let first_report = report.clone();
    let first_frame = Box::new(move || {
        times.first_frame_ms = Some(now_ms() - chosen_ms);
        let line = times.status(&name, cells, fields, &gpu);
        say(&line);
        tell(&first_report, &line);
    });
    let failed = Box::new(move |message: &str| {
        say(&format!("kuniumi-viewer: {message}"));
        tell(&report, &format!("error: {message}"));
    });
    let attributes = Window::default_attributes()
        .with_canvas(Some(canvas.clone()))
        .with_focusable(true)
        .with_prevent_default(true);
    app::spawn(scene, view, attributes, prepared, first_frame, failed)?;
    // Keys go to the canvas at once, without a click first.
    let _ = canvas.focus();
    Ok(status)
}
