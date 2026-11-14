//! The window: winit's event loop, the mouse and keyboard, and the
//! benchmark that spins the globe and times the frames.

use std::sync::Arc;
use std::time::Instant;

use kuniumi_math::units::RADIANS_PER_DEGREE;
use winit::application::ApplicationHandler;
use winit::dpi::{PhysicalPosition, PhysicalSize};
use winit::event::{ElementState, MouseButton, MouseScrollDelta, WindowEvent};
use winit::event_loop::{ActiveEventLoop, ControlFlow, EventLoop};
use winit::keyboard::{Key, NamedKey};
use winit::window::{Window, WindowId};

use crate::Scene;
use crate::gpu::{Renderer, block_on, request_device};
use crate::view::{View, relief_down, relief_up};

/// How far the benchmark turns the globe per frame (rad): 0.5°.
pub const BENCHMARK_STEP_RAD: f64 = 0.5 * RADIANS_PER_DEGREE;

/// Pixels of a touchpad scroll that count as one wheel step.
const PIXELS_PER_STEP: f64 = 40.0;

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
pub fn run(
    scene: Scene,
    view: View,
    size_px: (u32, u32),
    benchmark: Option<u32>,
) -> Result<(), String> {
    let event_loop = EventLoop::new().map_err(|e| format!("cannot open a window: {e}"))?;
    let mut app = App {
        scene,
        start: view,
        view,
        size_px,
        benchmark: benchmark.map(|frames| Benchmark {
            frames,
            presented: Vec::with_capacity(frames as usize + 1),
        }),
        gfx: None,
        error: None,
        dragging: false,
        cursor: None,
    };
    event_loop
        .run_app(&mut app)
        .map_err(|e| format!("the window failed: {e}"))?;
    app.error.map_or(Ok(()), Err)
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
    surface: wgpu::Surface<'static>,
    config: wgpu::SurfaceConfiguration,
    /// The sRGB view format of the surface's textures.
    view_format: wgpu::TextureFormat,
    renderer: Renderer,
}

struct App {
    scene: Scene,
    /// The view at start, restored by R.
    start: View,
    view: View,
    size_px: (u32, u32),
    benchmark: Option<Benchmark>,
    gfx: Option<Gfx>,
    error: Option<String>,
    dragging: bool,
    cursor: Option<PhysicalPosition<f64>>,
}

impl App {
    fn open(&self, event_loop: &ActiveEventLoop) -> Result<Gfx, String> {
        let attributes = Window::default_attributes()
            .with_title(self.title())
            .with_inner_size(PhysicalSize::new(self.size_px.0, self.size_px.1));
        let window = Arc::new(
            event_loop
                .create_window(attributes)
                .map_err(|e| format!("cannot open a window: {e}"))?,
        );
        let instance =
            wgpu::Instance::new(wgpu::InstanceDescriptor::new_with_display_handle_from_env(
                Box::new(event_loop.owned_display_handle()),
            ));
        let surface = instance
            .create_surface(Arc::clone(&window))
            .map_err(|e| format!("cannot draw into the window: {e}"))?;
        let adapter = block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: Some(&surface),
            ..Default::default()
        }))
        .map_err(|e| format!("no GPU adapter for the window: {e}"))?;
        let info = adapter.get_info();
        println!(
            "GPU: {} ({:?}, {:?})",
            info.name, info.backend, info.device_type
        );
        let (device, queue) = request_device(&adapter, &self.scene.mesh)?;
        let size = window.inner_size();
        let mut config = surface
            .get_default_config(&adapter, size.width.max(1), size.height.max(1))
            .ok_or("the GPU cannot draw into this window")?;
        let formats = surface.get_capabilities(&adapter).formats;
        config.format = formats
            .iter()
            .copied()
            .find(|f| f.is_srgb())
            .unwrap_or(config.format);
        let view_format = config.format.add_srgb_suffix();
        if view_format != config.format {
            config.view_formats = vec![view_format];
        }
        config.present_mode = wgpu::PresentMode::Fifo;
        surface.configure(&device, &config);
        let colours = self.scene.vertex_colours();
        let renderer = Renderer::new(
            &adapter,
            device,
            queue,
            view_format,
            &self.scene.mesh,
            &colours,
        );
        Ok(Gfx {
            window,
            surface,
            config,
            view_format,
            renderer,
        })
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
            self.view.relief,
            if self.view.sea { "on" } else { "off" }
        )
    }

    /// Shows the layer now selected: new colours, the legend printed.
    fn show_layer(&self) {
        if let Some(layer) = self.scene.layers.current() {
            println!(
                "layer {}/{}: {}",
                self.scene.layers.current_index() + 1,
                self.scene.layers.len(),
                layer.legend()
            );
        }
        if let Some(gfx) = &self.gfx {
            gfx.renderer.set_colours(&self.scene.vertex_colours());
        }
    }

    fn key(&mut self, key: &Key, event_loop: &ActiveEventLoop) {
        match press(key, &mut self.scene, &mut self.view, &self.start) {
            Action::Quit => event_loop.exit(),
            Action::Layer => {
                self.show_layer();
                self.changed();
            }
            Action::Redraw => self.changed(),
            Action::Nothing => {}
        }
    }

    /// Something visible changed: new title, new frame.
    fn changed(&self) {
        if let Some(gfx) = &self.gfx {
            gfx.window.set_title(&self.title());
            gfx.window.request_redraw();
        }
    }

    fn redraw(&mut self, event_loop: &ActiveEventLoop) {
        let Some(gfx) = &mut self.gfx else { return };
        let texture = match gfx.surface.get_current_texture() {
            wgpu::CurrentSurfaceTexture::Success(texture)
            | wgpu::CurrentSurfaceTexture::Suboptimal(texture) => texture,
            wgpu::CurrentSurfaceTexture::Timeout | wgpu::CurrentSurfaceTexture::Occluded => return,
            wgpu::CurrentSurfaceTexture::Outdated | wgpu::CurrentSurfaceTexture::Lost => {
                gfx.surface.configure(gfx.renderer.device(), &gfx.config);
                gfx.window.request_redraw();
                return;
            }
            wgpu::CurrentSurfaceTexture::Validation => {
                self.error = Some("the window's surface failed validation".into());
                event_loop.exit();
                return;
            }
        };
        let view = texture.texture.create_view(&wgpu::TextureViewDescriptor {
            format: Some(gfx.view_format),
            ..Default::default()
        });
        let (width, height) = (gfx.config.width, gfx.config.height);
        let frame = self
            .view
            .frame(f64::from(width) / f64::from(height), &self.scene.relief);
        gfx.renderer.draw(&view, width, height, &frame);
        gfx.window.pre_present_notify();
        gfx.renderer.present(texture);
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
    /// Not a key of the viewer.
    Nothing,
}

/// Applies the key `key` to the layers of `scene` and to `view` (`start` is
/// the view that R restores): Left/Right previous/next layer, Up/Down relief
/// steps, S sea on/off, R reset view, Esc or Q quit.
fn press(key: &Key, scene: &mut Scene, view: &mut View, start: &View) -> Action {
    let top = |view: &View| view.top(&scene.relief);
    match key {
        Key::Named(NamedKey::ArrowRight) => scene.layers.next(),
        Key::Named(NamedKey::ArrowLeft) => scene.layers.previous(),
        Key::Named(NamedKey::ArrowUp | NamedKey::ArrowDown) => {
            view.relief = if *key == Key::Named(NamedKey::ArrowUp) {
                relief_up(view.relief)
            } else {
                relief_down(view.relief)
            };
            view.camera.clamp(top(view));
            println!("relief x{}", view.relief);
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
                self.gfx = Some(gfx);
                self.show_layer();
                if self.benchmark.is_some() {
                    event_loop.set_control_flow(ControlFlow::Poll);
                }
                self.changed();
            }
            Err(e) => {
                self.error = Some(e);
                event_loop.exit();
            }
        }
    }

    fn window_event(&mut self, event_loop: &ActiveEventLoop, _: WindowId, event: WindowEvent) {
        match event {
            WindowEvent::CloseRequested => event_loop.exit(),
            WindowEvent::Resized(size) => {
                if let Some(gfx) = &mut self.gfx {
                    gfx.config.width = size.width.max(1);
                    gfx.config.height = size.height.max(1);
                    gfx.surface.configure(gfx.renderer.device(), &gfx.config);
                    gfx.window.request_redraw();
                }
            }
            WindowEvent::RedrawRequested => self.redraw(event_loop),
            WindowEvent::KeyboardInput { event, .. } if event.state == ElementState::Pressed => {
                self.key(&event.logical_key, event_loop);
            }
            WindowEvent::MouseInput {
                state,
                button: MouseButton::Left,
                ..
            } => self.dragging = state == ElementState::Pressed,
            WindowEvent::CursorMoved { position, .. } => {
                if let (true, Some(last), Some(gfx)) = (self.dragging, self.cursor, &self.gfx) {
                    let height = f64::from(gfx.config.height);
                    self.view
                        .camera
                        .drag(position.x - last.x, position.y - last.y, height);
                    self.changed();
                }
                self.cursor = Some(position);
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
    use kuniumi_world::{PlanetParams, Settings};

    use super::*;
    use crate::camera::OrbitCamera;

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
