//! The window: winit's event loop, the mouse and keyboard, the cell under
//! the cursor (highlighted, shown in the panel, printed on a click), and the
//! benchmark that spins the globe and times the frames.
//!
//! In a browser (ticket MV-03) the window is the page's canvas, whose GPU
//! side the page prepares first ([`prepare`], asynchronously); the event
//! loop is then spawned (`spawn`) rather than run, Q and Esc do nothing,
//! the title goes to the page's title and the printed lines to the console
//! ([`crate::say`]).

use std::sync::Arc;
use std::time::Instant;

use kuniumi_math::units::RADIANS_PER_DEGREE;
use winit::application::ApplicationHandler;
use winit::dpi::PhysicalPosition;
#[cfg(not(target_arch = "wasm32"))]
use winit::dpi::PhysicalSize;
use winit::event::{ElementState, MouseButton, MouseScrollDelta, WindowEvent};
use winit::event_loop::{ActiveEventLoop, ControlFlow, EventLoop};
use winit::keyboard::{Key, NamedKey};
use winit::window::{Window, WindowAttributes, WindowId};

#[cfg(not(target_arch = "wasm32"))]
use crate::gpu::block_on;
use crate::gpu::{Renderer, request_device, surface_config};
use crate::inspect::{NO_CELL, cell_text};
use crate::layers::highlight;
use crate::panel::{self, text_scale};
use crate::view::{View, relief_down, relief_line, relief_text, relief_up};
use crate::{Scene, say};

/// How far the benchmark turns the globe per frame (rad): 0.5°.
pub const BENCHMARK_STEP_RAD: f64 = 0.5 * RADIANS_PER_DEGREE;

/// Pixels of a touchpad scroll that count as one wheel step.
const PIXELS_PER_STEP: f64 = 40.0;

/// The largest distance between pressing and releasing the left button
/// (px) for a click, which prints the cell under the cursor; more is a drag.
pub const CLICK_PX: f64 = 3.0;

/// Whether the left button pressed at `press` and released at `release`
/// (pixel positions) is a click: at most [`CLICK_PX`] apart.
pub fn is_click(press: (f64, f64), release: (f64, f64)) -> bool {
    let (dx, dy) = (release.0 - press.0, release.1 - press.1);
    dx * dx + dy * dy <= CLICK_PX * CLICK_PX
}

/// Opens the window and runs until it is closed; with `benchmark`, spins the
/// globe for that many frames, prints the frame times and closes.
///
/// The window presents with vsync (`Fifo`). On Wayland, winit also paces
/// redraws by the compositor's frame callbacks, so no window exceeds the
/// display's refresh rate; `--benchmark --offscreen` times the drawing
/// itself.
///
/// # Errors
///
/// If no window or GPU can be opened.
#[cfg(not(target_arch = "wasm32"))]
pub fn run(
    scene: Scene,
    view: View,
    size_px: (u32, u32),
    benchmark: Option<u32>,
) -> Result<(), String> {
    let event_loop = EventLoop::new().map_err(|e| format!("cannot open a window: {e}"))?;
    let attributes =
        Window::default_attributes().with_inner_size(PhysicalSize::new(size_px.0, size_px.1));
    let mut app = App::new(scene, view, attributes, None);
    app.benchmark = benchmark.map(|frames| Benchmark {
        frames,
        presented: Vec::with_capacity(frames as usize + 1),
    });
    event_loop
        .run_app(&mut app)
        .map_err(|e| format!("the window failed: {e}"))?;
    app.error.map_or(Ok(()), Err)
}

/// Starts the event loop of the page's window, made with `attributes` (on
/// the canvas that `prepared` draws into), and returns at once: the
/// browser runs the loop. `first_frame` is called after the first frame is
/// presented; `failed` with the message of an error that stops the window.
///
/// # Errors
///
/// If winit's event loop cannot be made (a second one in the same page).
#[cfg(target_arch = "wasm32")]
pub fn spawn(
    scene: Scene,
    view: View,
    attributes: WindowAttributes,
    prepared: Prepared,
    first_frame: Box<dyn FnOnce()>,
    failed: Report,
) -> Result<(), String> {
    use winit::platform::web::EventLoopExtWebSys;
    let event_loop = EventLoop::new().map_err(|e| format!("cannot start the viewer: {e}"))?;
    let mut app = App::new(scene, view, attributes, Some(prepared));
    app.first_frame = Some(first_frame);
    app.failed = Some(failed);
    event_loop.spawn_app(app);
    Ok(())
}

/// A function told the message of an error that stops the window.
type Report = Box<dyn Fn(&str)>;

/// The GPU's side of a window: its surface, configured, and the globe's
/// renderer. In a browser it is prepared before the window exists, as it
/// needs to wait for the browser ([`prepare`]).
pub struct Prepared {
    surface: wgpu::Surface<'static>,
    config: wgpu::SurfaceConfiguration,
    /// The sRGB view format of the surface's textures.
    view_format: wgpu::TextureFormat,
    /// The largest side of the surface (px), the adapter's.
    max_side_px: u32,
    renderer: Renderer,
}

/// The GPU's side of a window for `scene`'s globe: the device of `adapter`,
/// `surface` configured for `width_px` × `height_px` ([`surface_config`])
/// and the renderer with the colours of the layer shown.
///
/// # Errors
///
/// If the device cannot be created or cannot draw into the surface.
pub async fn prepare(
    adapter: &wgpu::Adapter,
    surface: wgpu::Surface<'static>,
    width_px: u32,
    height_px: u32,
    scene: &Scene,
) -> Result<Prepared, String> {
    let (device, queue) = request_device(adapter, &scene.mesh).await?;
    let (config, view_format) = surface_config(&surface, adapter, width_px, height_px)?;
    surface.configure(&device, &config);
    let colours = scene.vertex_colours();
    let renderer = Renderer::new(adapter, device, queue, view_format, &scene.mesh, &colours);
    Ok(Prepared {
        surface,
        config,
        view_format,
        max_side_px: adapter.limits().max_texture_dimension_2d.max(1),
        renderer,
    })
}

/// The frame times of a benchmark.
struct Benchmark {
    frames: u32,
    /// When each frame was presented.
    presented: Vec<Instant>,
}

/// Summary of frame times.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameStats {
    /// Mean time between frames (ms).
    pub mean_ms: f64,
    /// Median (ms).
    pub median_ms: f64,
    /// 95th percentile, nearest rank (ms).
    pub p95_ms: f64,
    /// Longest (ms).
    pub max_ms: f64,
}

impl FrameStats {
    /// The statistics of the times between frames `intervals_ms` (at least
    /// one).
    pub fn of(intervals_ms: &[f64]) -> Self {
        let mut sorted = intervals_ms.to_vec();
        sorted.sort_by(f64::total_cmp);
        let n = sorted.len();
        // Nearest rank: the smallest value with at least p·n values at or
        // below it.
        let rank = |p: f64| sorted[((p * n as f64).ceil() as usize).clamp(1, n) - 1];
        FrameStats {
            mean_ms: sorted.iter().sum::<f64>() / n as f64,
            median_ms: rank(0.5),
            p95_ms: rank(0.95),
            max_ms: sorted[n - 1],
        }
    }

    /// Prints the statistics of `n` frames of `triangles` triangles drawn as
    /// `what` says.
    pub fn print(&self, what: &str, triangles: usize) {
        println!(
            "benchmark, {triangles} triangles, {what}: mean {:.2} ms ({:.1} fps), median {:.2} ms, \
             p95 {:.2} ms, max {:.2} ms",
            self.mean_ms,
            1e3 / self.mean_ms,
            self.median_ms,
            self.p95_ms,
            self.max_ms
        );
    }
}

/// The window and its GPU state.
struct Gfx {
    window: Arc<Window>,
    gpu: Prepared,
}

struct App {
    scene: Scene,
    /// The view at start, restored by R.
    start: View,
    view: View,
    /// The window's attributes, until it is made.
    attributes: Option<WindowAttributes>,
    /// The GPU's side of the window, if prepared before it (in a browser).
    prepared: Option<Prepared>,
    benchmark: Option<Benchmark>,
    gfx: Option<Gfx>,
    /// The error that stopped the window, returned by `run`.
    error: Option<String>,
    /// Called after the first frame is presented (in a browser).
    first_frame: Option<Box<dyn FnOnce()>>,
    /// Called with an error that stops the window (in a browser, where
    /// nothing returns it).
    failed: Option<Report>,
    dragging: bool,
    /// The cursor's position in the window, if it is inside.
    cursor: Option<PhysicalPosition<f64>>,
    /// Where the left button was last pressed.
    pressed_at: Option<PhysicalPosition<f64>>,
    /// The cell under the cursor, highlighted.
    picked: Option<u32>,
    /// Whether the panel is shown (P).
    panel: bool,
    /// The panel's text scale, from the window's scale factor.
    text_scale: u32,
}

impl App {
    fn new(
        scene: Scene,
        view: View,
        attributes: WindowAttributes,
        prepared: Option<Prepared>,
    ) -> Self {
        App {
            scene,
            start: view,
            view,
            attributes: Some(attributes),
            prepared,
            benchmark: None,
            gfx: None,
            error: None,
            first_frame: None,
            failed: None,
            dragging: false,
            cursor: None,
            pressed_at: None,
            picked: None,
            panel: true,
            text_scale: 1,
        }
    }

    fn open(&mut self, event_loop: &ActiveEventLoop) -> Result<Gfx, String> {
        let attributes = self
            .attributes
            .take()
            .ok_or("the window was already made")?
            .with_title(self.title());
        let window = Arc::new(
            event_loop
                .create_window(attributes)
                .map_err(|e| format!("cannot open a window: {e}"))?,
        );
        let gpu = match self.prepared.take() {
            Some(gpu) => gpu,
            #[cfg(not(target_arch = "wasm32"))]
            None => self.prepare_window(event_loop, &window)?,
            #[cfg(target_arch = "wasm32")]
            None => return Err("no GPU was prepared for the canvas".into()),
        };
        Ok(Gfx { window, gpu })
    }

    /// The GPU's side of the native `window`.
    #[cfg(not(target_arch = "wasm32"))]
    fn prepare_window(
        &self,
        event_loop: &ActiveEventLoop,
        window: &Arc<Window>,
    ) -> Result<Prepared, String> {
        let instance =
            wgpu::Instance::new(wgpu::InstanceDescriptor::new_with_display_handle_from_env(
                Box::new(event_loop.owned_display_handle()),
            ));
        let surface = instance
            .create_surface(Arc::clone(window))
            .map_err(|e| format!("cannot draw into the window: {e}"))?;
        let adapter = block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: Some(&surface),
            ..Default::default()
        }))
        .map_err(|e| format!("no GPU adapter for the window: {e}"))?;
        let info = adapter.get_info();
        say(&format!(
            "GPU: {} ({:?}, {:?})",
            info.name, info.backend, info.device_type
        ));
        let size = window.inner_size();
        block_on(prepare(
            &adapter,
            surface,
            size.width,
            size.height,
            &self.scene,
        ))
    }

    /// Stops the window with the error `message`.
    fn fail(&mut self, event_loop: &ActiveEventLoop, message: String) {
        if let Some(failed) = &self.failed {
            failed(&message);
        }
        self.error = Some(message);
        event_loop.exit();
    }

    fn title(&self) -> String {
        let layer = match self.scene.layers.current() {
            Some(layer) => format!(
                "{}/{} {} [{}]",
                self.scene.layers.current_index() + 1,
                self.scene.layers.len(),
                layer.name,
                if layer.unit.is_empty() {
                    "-"
                } else {
                    &layer.unit
                }
            ),
            None => "no fields".into(),
        };
        format!(
            "kuniumi-viewer - {} - {layer} - relief x{} - sea {}",
            self.scene.dir.display(),
            relief_text(self.view.relief),
            if self.view.sea { "on" } else { "off" }
        )
    }

    /// Shows the layer now selected: new colours, the legend printed, the
    /// panel's legend.
    fn show_layer(&mut self) {
        if let Some(layer) = self.scene.layers.current() {
            say(&format!(
                "layer {}/{}: {}",
                self.scene.layers.current_index() + 1,
                self.scene.layers.len(),
                layer.legend()
            ));
        }
        if let Some(gfx) = &self.gfx {
            gfx.gpu.renderer.set_colours(&self.scene.vertex_colours());
        }
        if let Some(cell) = self.picked {
            self.paint(cell, true);
        }
        self.update_panel();
    }

    /// Lays out the panel anew (none while it is hidden) for the cell under
    /// the cursor, the layer, the window's size and its text scale, and
    /// asks for a new frame.
    fn update_panel(&mut self) {
        let Some(gfx) = &mut self.gfx else { return };
        let quads = if self.panel {
            let lines = panel::lines(&self.scene.world, &self.scene.layers, self.picked);
            let viewport = (gfx.gpu.config.width, gfx.gpu.config.height);
            panel::layout(
                &lines,
                self.scene.layers.current(),
                self.text_scale,
                viewport,
            )
        } else {
            Vec::new()
        };
        gfx.gpu.renderer.set_overlay(&quads);
        gfx.window.request_redraw();
    }

    /// Paints `cell` in the layer's colour, highlighted or not.
    fn paint(&self, cell: u32, highlighted: bool) {
        if let Some(gfx) = &self.gfx {
            let colour = self.scene.layers.cell_colour(&self.scene.world, cell);
            let colour = if highlighted {
                highlight(colour)
            } else {
                colour
            };
            gfx.gpu
                .renderer
                .set_cell_colour(self.scene.mesh.cell_vertices(cell), colour);
        }
    }

    /// Picks the cell under the cursor anew, after the cursor moved or the
    /// view changed, and moves the highlight and the panel to it (none while
    /// the benchmark runs).
    fn repick(&mut self) {
        let Some(gfx) = &self.gfx else { return };
        let (width, height) = (gfx.gpu.config.width, gfx.gpu.config.height);
        let cell = match (self.cursor, &self.benchmark) {
            (Some(cursor), None) => self
                .scene
                .pick(&self.view, width, height, cursor.x, cursor.y)
                .map(|hit| hit.cell),
            _ => None,
        };
        if cell != self.picked {
            if let Some(old) = self.picked {
                self.paint(old, false);
            }
            if let Some(new) = cell {
                self.paint(new, true);
            }
            self.picked = cell;
            self.update_panel();
        }
    }

    /// Prints the cell under the cursor to the terminal (a click).
    fn print_picked(&self) {
        match self.picked {
            Some(cell) => {
                for line in cell_text(&self.scene.world, &self.scene.layers, cell) {
                    say(&line);
                }
            }
            None => say(NO_CELL),
        }
    }

    fn key(&mut self, key: &Key, event_loop: &ActiveEventLoop) {
        match press(key, &mut self.scene, &mut self.view, &self.start) {
            // A page is closed by the browser, not by a key.
            Action::Quit if cfg!(target_arch = "wasm32") => {}
            Action::Quit => event_loop.exit(),
            Action::Layer => {
                self.show_layer();
                self.changed();
            }
            Action::Redraw => {
                self.changed();
                self.repick();
            }
            Action::Panel => {
                self.panel = !self.panel;
                self.update_panel();
            }
            Action::Nothing => {}
        }
    }

    /// Something visible changed: new title, new frame.
    fn changed(&self) {
        if let Some(gfx) = &self.gfx {
            let title = self.title();
            gfx.window.set_title(&title);
            // winit puts a canvas's title into its `alt` attribute only.
            #[cfg(target_arch = "wasm32")]
            crate::web::set_title(&title);
            gfx.window.request_redraw();
        }
    }

    fn redraw(&mut self, event_loop: &ActiveEventLoop) {
        let Some(gfx) = &mut self.gfx else { return };
        let texture = match gfx.gpu.surface.get_current_texture() {
            wgpu::CurrentSurfaceTexture::Success(texture)
            | wgpu::CurrentSurfaceTexture::Suboptimal(texture) => texture,
            wgpu::CurrentSurfaceTexture::Timeout | wgpu::CurrentSurfaceTexture::Occluded => return,
            wgpu::CurrentSurfaceTexture::Outdated | wgpu::CurrentSurfaceTexture::Lost => {
                gfx.gpu
                    .surface
                    .configure(gfx.gpu.renderer.device(), &gfx.gpu.config);
                gfx.window.request_redraw();
                return;
            }
            wgpu::CurrentSurfaceTexture::Validation => {
                self.fail(event_loop, "the window's surface failed validation".into());
                return;
            }
        };
        let view = texture.texture.create_view(&wgpu::TextureViewDescriptor {
            format: Some(gfx.gpu.view_format),
            ..Default::default()
        });
        let (width, height) = (gfx.gpu.config.width, gfx.gpu.config.height);
        let frame = self
            .view
            .frame(f64::from(width) / f64::from(height), &self.scene.relief);
        gfx.gpu.renderer.draw(&view, width, height, &frame);
        gfx.window.pre_present_notify();
        gfx.gpu.renderer.present(texture);
        if let Some(first_frame) = self.first_frame.take() {
            first_frame();
        }
        if let Some(benchmark) = &mut self.benchmark {
            benchmark.presented.push(Instant::now());
            self.view.camera.lon_rad += BENCHMARK_STEP_RAD;
            if benchmark.presented.len() > benchmark.frames as usize {
                let intervals: Vec<f64> = benchmark
                    .presented
                    .windows(2)
                    .map(|pair| (pair[1] - pair[0]).as_secs_f64() * 1e3)
                    .collect();
                let what = format!("{width}x{height} px in the window, between presents");
                FrameStats::of(&intervals).print(&what, self.scene.mesh.indices.len() / 3);
                event_loop.exit();
            }
        }
    }
}

/// What a key press asks of the window.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Action {
    /// Close the window.
    Quit,
    /// Another layer: new colours and a new frame.
    Layer,
    /// The view changed: a new frame.
    Redraw,
    /// Show or hide the panel.
    Panel,
    /// Not a key of the viewer.
    Nothing,
}

/// Applies the key `key` to the layers of `scene` and to `view` (`start` is
/// the view that R restores): Left/Right previous/next layer, Up/Down relief
/// steps, S sea on/off, R reset view, P panel on/off (left to the window),
/// Esc or Q quit.
fn press(key: &Key, scene: &mut Scene, view: &mut View, start: &View) -> Action {
    let top = |view: &View| view.top(&scene.relief);
    match key {
        Key::Named(NamedKey::ArrowRight) => scene.layers.next(),
        Key::Named(NamedKey::ArrowLeft) => scene.layers.previous(),
        Key::Named(NamedKey::ArrowUp | NamedKey::ArrowDown) => {
            view.relief = if *key == Key::Named(NamedKey::ArrowUp) {
                relief_up(view.relief, scene.relief.max_relief())
            } else {
                relief_down(view.relief)
            };
            view.camera.clamp(top(view));
            say(&relief_line(view.relief, scene.relief.max_relief()));
            return Action::Redraw;
        }
        Key::Named(NamedKey::Escape) => return Action::Quit,
        Key::Character(c) => {
            match c.as_str() {
                "s" | "S" => view.sea = !view.sea,
                "r" | "R" => {
                    view.camera = start.camera;
                    view.camera.clamp(top(view));
                }
                "q" | "Q" => return Action::Quit,
                "p" | "P" => return Action::Panel,
                _ => return Action::Nothing,
            }
            return Action::Redraw;
        }
        _ => return Action::Nothing,
    }
    Action::Layer
}

impl ApplicationHandler for App {
    fn resumed(&mut self, event_loop: &ActiveEventLoop) {
        if self.gfx.is_some() {
            return;
        }
        match self.open(event_loop) {
            Ok(gfx) => {
                self.text_scale = text_scale(gfx.window.scale_factor());
                self.gfx = Some(gfx);
                self.show_layer();
                if self.benchmark.is_some() {
                    event_loop.set_control_flow(ControlFlow::Poll);
                }
                self.changed();
            }
            Err(e) => self.fail(event_loop, e),
        }
    }

    fn window_event(&mut self, event_loop: &ActiveEventLoop, _: WindowId, event: WindowEvent) {
        match event {
            WindowEvent::CloseRequested => event_loop.exit(),
            WindowEvent::Resized(size) => {
                if let Some(gfx) = &mut self.gfx {
                    let side = gfx.gpu.max_side_px;
                    gfx.gpu.config.width = size.width.clamp(1, side);
                    gfx.gpu.config.height = size.height.clamp(1, side);
                    gfx.gpu
                        .surface
                        .configure(gfx.gpu.renderer.device(), &gfx.gpu.config);
                    gfx.window.request_redraw();
                }
                self.repick();
                self.update_panel();
            }
            WindowEvent::ScaleFactorChanged { scale_factor, .. } => {
                self.text_scale = text_scale(scale_factor);
                self.update_panel();
            }
            WindowEvent::RedrawRequested => self.redraw(event_loop),
            WindowEvent::KeyboardInput { event, .. } if event.state == ElementState::Pressed => {
                self.key(&event.logical_key, event_loop);
            }
            WindowEvent::MouseInput {
                state,
                button: MouseButton::Left,
                ..
            } => {
                self.dragging = state == ElementState::Pressed;
                if self.dragging {
                    self.pressed_at = self.cursor;
                } else if let (Some(press), Some(release)) = (self.pressed_at.take(), self.cursor) {
                    if is_click((press.x, press.y), (release.x, release.y)) {
                        self.print_picked();
                    }
                }
            }
            WindowEvent::CursorMoved { position, .. } => {
                if let (true, Some(last), Some(gfx)) = (self.dragging, self.cursor, &self.gfx) {
                    let height = f64::from(gfx.gpu.config.height);
                    self.view
                        .camera
                        .drag(position.x - last.x, position.y - last.y, height);
                    self.changed();
                }
                self.cursor = Some(position);
                self.repick();
            }
            WindowEvent::CursorLeft { .. } => {
                self.cursor = None;
                self.repick();
            }
            WindowEvent::MouseWheel { delta, .. } => {
                let steps = match delta {
                    MouseScrollDelta::LineDelta(_, y) => f64::from(y),
                    MouseScrollDelta::PixelDelta(p) => p.y / PIXELS_PER_STEP,
                };
                self.view
                    .camera
                    .zoom(steps, self.view.top(&self.scene.relief));
                self.changed();
                self.repick();
            }
            _ => {}
        }
    }

    fn about_to_wait(&mut self, _: &ActiveEventLoop) {
        if let (Some(_), Some(gfx)) = (&self.benchmark, &self.gfx) {
            gfx.window.request_redraw();
        }
    }
}

#[cfg(test)]
mod tests {
    use kuniumi_sim::stages::crust::SURFACE_HEIGHT;
    use kuniumi_sim::{StageId, pipeline};
    use kuniumi_world::{PlanetParams, Settings};

    use super::*;
    use crate::args::Parsed;
    use crate::camera::OrbitCamera;
    use crate::start_view;
    use crate::view::{MAX_RELIEF, MIN_DRAWN_RADIUS};

    #[test]
    fn keys() {
        let settings = Settings {
            cells: 1_000,
            allow_partial: true,
            ..Settings::default()
        };
        let world =
            kuniumi_sim::generate(&PlanetParams::earth(), &settings, 3, &mut |_| {}).unwrap();
        let mut scene = Scene::new("w".into(), world, Some("plate_id")).unwrap();
        let start = View {
            camera: OrbitCamera {
                lat_rad: 0.2,
                lon_rad: 0.0,
                distance: 3.0,
            },
            relief: 20.0,
            sea: true,
        };
        let mut view = start;
        let mut hit = |key: Key, view: &mut View| press(&key, &mut scene, view, &start);
        let named = Key::Named;
        let character = |c: &str| Key::Character(c.into());
        assert_eq!(hit(named(NamedKey::ArrowRight), &mut view), Action::Layer);
        assert_eq!(hit(named(NamedKey::ArrowLeft), &mut view), Action::Layer);
        assert_eq!(hit(named(NamedKey::ArrowLeft), &mut view), Action::Layer);
        assert_eq!(hit(named(NamedKey::ArrowUp), &mut view), Action::Redraw);
        assert_eq!(view.relief.to_bits(), 50.0_f64.to_bits());
        assert_eq!(hit(named(NamedKey::ArrowDown), &mut view), Action::Redraw);
        assert_eq!(hit(named(NamedKey::ArrowDown), &mut view), Action::Redraw);
        assert_eq!(view.relief.to_bits(), 10.0_f64.to_bits());
        assert_eq!(hit(character("s"), &mut view), Action::Redraw);
        assert!(!view.sea);
        assert_eq!(hit(character("S"), &mut view), Action::Redraw);
        assert!(view.sea);
        view.camera.drag(300.0, -200.0, 800.0);
        view.camera.zoom(5.0, 1.0);
        assert_ne!(view.camera, start.camera);
        assert_eq!(hit(character("r"), &mut view), Action::Redraw);
        assert_eq!(view.camera, start.camera);
        let before = view;
        assert_eq!(hit(character("p"), &mut view), Action::Panel);
        assert_eq!(hit(character("P"), &mut view), Action::Panel);
        assert_eq!(view, before, "the panel is not part of the view");
        assert_eq!(hit(character("x"), &mut view), Action::Nothing);
        assert_eq!(hit(named(NamedKey::Tab), &mut view), Action::Nothing);
        assert_eq!(hit(character("q"), &mut view), Action::Quit);
        assert_eq!(hit(named(NamedKey::Escape), &mut view), Action::Quit);
        // Right, Left, Left: one layer before plate_id.
        let plate_id = scene
            .world
            .fields
            .iter()
            .position(|f| f.name == "plate_id")
            .unwrap();
        assert_eq!(scene.layers.current_index(), plate_id - 1);
        // A relief too high for the camera pushes it out.
        let mut close = View {
            camera: OrbitCamera {
                distance: 1.0 + 1e-3,
                ..start.camera
            },
            relief: 50.0,
            ..start
        };
        assert_eq!(
            press(&named(NamedKey::ArrowUp), &mut scene, &mut close, &start),
            Action::Redraw
        );
        assert_eq!(
            close.camera.distance.to_bits(),
            (close.top(&scene.relief) + crate::camera::MIN_CLEARANCE).to_bits()
        );
    }

    /// The planet of review R3-1 of MV-02: the smallest radius of D-12 with
    /// Earth's bulk density (g = 9.80665 m/s² × 1,500/6,371 = 2.3089 m/s²)
    /// and, like every other parameter, Earth's ocean volume.
    fn small_planet() -> PlanetParams {
        PlanetParams {
            radius_m: 1_500_000.0,
            surface_gravity_m_s2: 2.3089,
            ..PlanetParams::earth()
        }
    }

    /// The scene of `params`, seed 1, 100,000 cells, through the Crust
    /// stage (its relief: the deepest point, which bounds the relief, is
    /// ocean floor that erosion leaves alone; the climate and 1,000 erosion
    /// steps at 100,000 cells would take many minutes in a debug build).
    fn scene_of(params: &PlanetParams) -> Scene {
        let settings = Settings {
            cells: 100_000,
            allow_partial: true,
            ..Settings::default()
        };
        let through_crust = |stage: StageId| {
            if stage <= StageId::Crust {
                pipeline::implementation(stage)
            } else {
                pipeline::Implementation::Missing
            }
        };
        let world = pipeline::run(params, &settings, 1, &through_crust, &mut |_| {}).unwrap();
        Scene::new("w".into(), world, None).unwrap()
    }

    /// Every relief exaggeration the viewer reaches on `scene`, without
    /// repeats: each start relief that `--relief` gives (and the default),
    /// through `start_view`, and each relief that ten presses of Up, then
    /// ten of Down, then ten of Up reach from it.
    fn reachable_reliefs(scene: &mut Scene) -> Vec<f64> {
        let mut lists: Vec<Vec<String>> = vec![vec!["w".into()]];
        let given = ["0", "0.5", "14", "14.79", "15", "37", "99.5"];
        let steps = crate::view::RELIEF_STEPS.map(|step| step.to_string());
        for relief in given.iter().map(|r| r.to_string()).chain(steps) {
            lists.push(vec!["w".into(), "--relief".into(), relief]);
        }
        let mut reached = Vec::new();
        for list in lists {
            let Ok(Parsed::Run(args)) = crate::args::parse(list.clone()) else {
                panic!("{list:?}");
            };
            let start = start_view(&args, &scene.relief);
            reached.push(start.relief);
            let mut view = start;
            for key in [NamedKey::ArrowUp, NamedKey::ArrowDown, NamedKey::ArrowUp] {
                for _ in 0..10 {
                    press(&Key::Named(key), scene, &mut view, &start);
                    reached.push(view.relief);
                }
            }
        }
        reached.sort_by(f64::total_cmp);
        reached.dedup_by(|a, b| a.to_bits() == b.to_bits());
        reached
    }

    /// The lowest radius at which a vertex of `scene` is drawn at `relief`,
    /// asserting that every vertex lies at [`MIN_DRAWN_RADIUS`] or above:
    /// in `f32` as the shader computes `1 + σ·h`, with separate roundings
    /// and with a fused multiply-add, and in `f64` from the `f32` values as
    /// the picker does, with the `σ` the frame hands to the GPU.
    fn lowest_drawn_radius(scene: &Scene, relief: f64) -> f64 {
        let view = View {
            camera: OrbitCamera {
                lat_rad: 0.0,
                lon_rad: 0.0,
                distance: 3.0,
            },
            relief,
            sea: true,
        };
        let sigma = view.frame(1.6, &scene.relief).relief_scale_per_m;
        let floor = MIN_DRAWN_RADIUS as f32;
        let mut lowest = f64::INFINITY;
        for (k, v) in scene.mesh.vertices.iter().enumerate() {
            let h = v[3];
            let separate = 1.0_f32 + sigma * h;
            let fused = sigma.mul_add(h, 1.0);
            let picker = 1.0 + f64::from(sigma) * f64::from(h);
            assert!(
                separate >= floor && fused >= floor && picker >= MIN_DRAWN_RADIUS,
                "×{relief}, vertex {k} at {h} m: {separate}, {fused}, {picker}"
            );
            if picker < lowest {
                lowest = picker;
            }
        }
        lowest
    }

    #[test]
    fn relief_keeps_the_deepest_point_of_a_small_world_above_the_floor() {
        let mut scene = scene_of(&small_planet());
        // The deepest cell below the sea level, from the world's fields.
        let world = &scene.world;
        let heights = &world.field(SURFACE_HEIGHT).unwrap().values;
        let depth_m = heights
            .iter()
            .map(|h| world.sea_level_m - h)
            .fold(0.0, f64::max);
        assert!(
            (50_000.0..51_000.0).contains(&depth_m),
            "50.7 km in review R3-1 of MV-02: {depth_m}"
        );
        // By hand: (1 − 0.5) × 1,500 km/50.7 km = 14.80.
        let cap = scene.relief.max_relief();
        let by_hand = 0.5 * 1_500_000.0 / depth_m;
        assert!((cap - by_hand).abs() < 1e-4, "{cap} {by_hand}");
        assert!(cap < MAX_RELIEF);
        let reached = reachable_reliefs(&mut scene);
        assert_eq!(reached.last().unwrap().to_bits(), cap.to_bits());
        for &relief in &reached {
            lowest_drawn_radius(&scene, relief);
        }
        // At the cap the deepest point lies just above half the radius.
        let lowest = lowest_drawn_radius(&scene, cap);
        assert!(lowest < MIN_DRAWN_RADIUS + 1e-5, "{lowest}");
        // Above the cap it would not (the bound is tight, not a margin).
        let view = View {
            camera: OrbitCamera {
                lat_rad: 0.0,
                lon_rad: 0.0,
                distance: 3.0,
            },
            relief: 15.0,
            sea: true,
        };
        let sigma = f64::from(view.frame(1.6, &scene.relief).relief_scale_per_m);
        assert!(1.0 - sigma * depth_m < MIN_DRAWN_RADIUS);
    }

    #[test]
    fn relief_on_earth_reaches_the_largest_step_above_the_floor() {
        let mut scene = scene_of(&PlanetParams::earth());
        assert_eq!(scene.relief.max_relief().to_bits(), MAX_RELIEF.to_bits());
        let reached = reachable_reliefs(&mut scene);
        assert_eq!(reached.last().unwrap().to_bits(), MAX_RELIEF.to_bits());
        for &relief in &reached {
            lowest_drawn_radius(&scene, relief);
        }
        // At ×100 Earth's deepest cell, some 7 km down, lies near 0.89.
        let lowest = lowest_drawn_radius(&scene, MAX_RELIEF);
        assert!((0.85..0.92).contains(&lowest), "{lowest}");
    }

    #[test]
    fn clicks_and_drags() {
        assert!(is_click((100.0, 200.0), (100.0, 200.0)));
        assert!(is_click((100.0, 200.0), (103.0, 200.0)));
        assert!(is_click((100.0, 200.0), (100.0, 197.0)));
        // 3-4-5: 1.8 and 2.4 px are 3 px apart.
        assert!(is_click((0.0, 0.0), (1.8, 2.4)));
        assert!(!is_click((100.0, 200.0), (103.1, 200.0)));
        assert!(!is_click((0.0, 0.0), (2.2, 2.2)));
    }

    #[test]
    fn frame_statistics_by_hand() {
        // 20 intervals: 1 … 20 ms. Mean 10.5; median the 10th, 10; the 95th
        // percentile the 19th, 19; maximum 20.
        let intervals: Vec<f64> = (1..=20).rev().map(f64::from).collect();
        let stats = FrameStats::of(&intervals);
        assert_eq!(
            stats,
            FrameStats {
                mean_ms: 10.5,
                median_ms: 10.0,
                p95_ms: 19.0,
                max_ms: 20.0
            }
        );
        let one = FrameStats::of(&[16.7]);
        assert_eq!(
            (one.mean_ms, one.median_ms, one.p95_ms, one.max_ms),
            (16.7, 16.7, 16.7, 16.7)
        );
    }
}
