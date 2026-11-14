//! `kuniumi-viewer`: a world directory as an interactive 3D globe (ticket
//! MV-02; D-19).
//!
//! `kuniumi-viewer <world dir>` reads the world (`kuniumi_io::world_dir`),
//! builds a triangle mesh of its Voronoi cells ([`mesh`]) and shows any of
//! its fields ([`layers`]) in the colours of the world's PNG maps, with the
//! relief exaggerated and a translucent sea-level surface ([`view`]), on the
//! GPU through `wgpu` ([`gpu`]) in a `winit` window ([`app`]). The mouse
//! turns and zooms the globe ([`camera`]), the keyboard changes the layer,
//! the relief and the sea. The cell under the cursor is found on the drawn
//! triangles ([`pick`]) and highlighted; a click prints its values
//! ([`inspect`]). `--screenshot` renders one frame into a PNG without a
//! window (with `--cursor`, the cell under a pixel picked, highlighted and
//! printed); `--benchmark` times the frames of a spinning globe. The usage
//! is in [`args::USAGE`].

// The GPU works in `f32`. D-19 and rule 1 of docs/conventions.md,
// "Determinism", let this crate (and kuniumi-gpu) use `f32` through this
// crate-level allowance: nothing the viewer computes feeds back into
// generation, and no other crate may do this.
#![allow(clippy::disallowed_types)]

mod app;
mod args;
mod camera;
mod gpu;
mod inspect;
mod layers;
mod mesh;
mod pick;
mod view;

use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Instant;

use kuniumi_io::Rgb8;
use kuniumi_sim::stages::crust::SURFACE_HEIGHT;
use kuniumi_world::World;

use crate::app::FrameStats;
use crate::args::{Args, Mode, Parsed, USAGE};
use crate::camera::{OrbitCamera, cursor_ndc};
use crate::gpu::Offscreen;
use crate::layers::Layers;
use crate::mesh::GlobeMesh;
use crate::pick::{Hit, Picker};
use crate::view::{Relief, View};

/// A world and what the viewer derives from it once.
pub struct Scene {
    /// The world directory.
    pub dir: PathBuf,
    /// The world.
    pub world: World,
    /// Its globe.
    pub mesh: GlobeMesh,
    /// Its layers, one per field.
    pub layers: Layers,
    /// Its relief as drawn.
    pub relief: Relief,
    /// Finds the cell under the cursor.
    pub picker: Picker,
}

impl Scene {
    /// The scene of `world` from the directory `dir`, showing the field
    /// `layer` if given.
    ///
    /// # Errors
    ///
    /// If the planet's radius is not finite and positive (a hand-edited
    /// `world.json`: `world_dir::read` does not check the parameters), or
    /// the world has no field `layer`.
    pub fn new(dir: PathBuf, world: World, layer: Option<&str>) -> Result<Self, String> {
        let radius_m = world.params.radius_m;
        if !(radius_m.is_finite() && radius_m > 0.0) {
            return Err(format!(
                "the planet's radius is {radius_m} m; the globe needs a finite positive radius"
            ));
        }
        let heights = world.field(SURFACE_HEIGHT).map(|field| &field.values);
        let above_sea: Vec<f64> = match heights {
            Some(h) => h.iter().map(|h| h - world.sea_level_m).collect(),
            None => vec![0.0; world.grid.len()],
        };
        let mesh = GlobeMesh::new(&world.grid, &above_sea);
        let mut layers = Layers::new(&world);
        if let Some(name) = layer {
            layers.select(name)?;
        }
        let relief = Relief {
            radius_m,
            max_height_m: mesh.max_height_m,
            present: heights.is_some(),
        };
        let picker = Picker::new(&mesh);
        Ok(Self {
            dir,
            world,
            mesh,
            layers,
            relief,
            picker,
        })
    }

    /// The colour of every vertex in the layer shown.
    pub fn vertex_colours(&self) -> Vec<[u8; 4]> {
        self.mesh.vertex_colours(&self.layers.colours(&self.world))
    }

    /// The cell under the point `(x_px, y_px)` (from the top-left corner) of
    /// a viewport `width_px` × `height_px` showing `view`, as drawn.
    pub fn pick(
        &mut self,
        view: &View,
        width_px: u32,
        height_px: u32,
        x_px: f64,
        y_px: f64,
    ) -> Option<Hit> {
        let aspect = f64::from(width_px) / f64::from(height_px);
        let ndc = cursor_ndc(x_px, y_px, width_px, height_px);
        let (origin, dir) = view.camera.ray(aspect, ndc);
        // The σ of the frame, as the GPU gets it.
        let sigma = view.relief_scale_per_m(&self.relief) as f32;
        self.picker
            .pick(&self.world.grid, &self.mesh, sigma, origin, dir)
    }
}

/// The length of each arm of [`crosshair`] (px), from the pixel to its end.
const CROSSHAIR_END_PX: i64 = 12;
/// The distance from the marked pixel at which the arms start (px).
const CROSSHAIR_START_PX: i64 = 3;

/// Marks the pixel `(x_px, y_px)` of the row-major image `pixels`
/// (`width_px` × `height_px`) with a crosshair: four white arms from
/// [`CROSSHAIR_START_PX`] to [`CROSSHAIR_END_PX`] pixels away, outlined in
/// black, clipped at the image's edges. The pixel and its eight neighbours
/// are left as drawn, so the cell under it stays visible.
pub fn crosshair(pixels: &mut [Rgb8], width_px: u32, height_px: u32, x_px: u32, y_px: u32) {
    const WHITE: Rgb8 = [255, 255, 255];
    const BLACK: Rgb8 = [0, 0, 0];
    let (w, h) = (i64::from(width_px), i64::from(height_px));
    let (x0, y0) = (i64::from(x_px), i64::from(y_px));
    let mut paint = |x: i64, y: i64, colour: Rgb8| {
        if (0..w).contains(&x) && (0..h).contains(&y) {
            pixels[(y * w + x) as usize] = colour;
        }
    };
    let arms = [(1, 0), (-1, 0), (0, 1), (0, -1)];
    // The outline first: both sides of each arm and a cap at either end.
    for (dx, dy) in arms {
        for r in CROSSHAIR_START_PX - 1..=CROSSHAIR_END_PX + 1 {
            let (x, y) = (x0 + r * dx, y0 + r * dy);
            paint(x + dy, y + dx, BLACK);
            paint(x - dy, y - dx, BLACK);
        }
        for r in [CROSSHAIR_START_PX - 1, CROSSHAIR_END_PX + 1] {
            paint(x0 + r * dx, y0 + r * dy, BLACK);
        }
    }
    for (dx, dy) in arms {
        for r in CROSSHAIR_START_PX..=CROSSHAIR_END_PX {
            paint(x0 + r * dx, y0 + r * dy, WHITE);
        }
    }
}

fn main() -> ExitCode {
    match args::parse(std::env::args().skip(1)) {
        Err(message) => {
            eprintln!("kuniumi-viewer: {message}\n\n{USAGE}");
            ExitCode::from(2)
        }
        Ok(Parsed::Help) => {
            println!("{USAGE}");
            ExitCode::SUCCESS
        }
        Ok(Parsed::Run(args)) => match run(&args) {
            Ok(()) => ExitCode::SUCCESS,
            Err(message) => {
                eprintln!("kuniumi-viewer: {message}");
                ExitCode::FAILURE
            }
        },
    }
}

fn run(args: &Args) -> Result<(), String> {
    let start = Instant::now();
    let world = kuniumi_io::world_dir::read(&args.world_dir).map_err(|e| e.to_string())?;
    let read_s = start.elapsed().as_secs_f64();
    let mut scene = Scene::new(args.world_dir.clone(), world, args.layer.as_deref())?;
    println!(
        "{}: {} cells, {} fields; read in {read_s:.2} s, mesh of {} triangles, layers and picker \
         in {:.2} s",
        args.world_dir.display(),
        scene.world.grid.len(),
        scene.world.fields.len(),
        scene.mesh.indices.len() / 3,
        start.elapsed().as_secs_f64() - read_s
    );
    let mut view = View {
        camera: OrbitCamera {
            lat_rad: args.lat_rad,
            lon_rad: args.lon_rad,
            distance: args.distance,
        },
        relief: args.relief,
        sea: args.sea,
    };
    view.camera.clamp(view.top(&scene.relief));
    let size = (args.width_px, args.height_px);
    let aspect = f64::from(size.0) / f64::from(size.1);
    match &args.mode {
        Mode::Screenshot(path) => {
            if let Some(layer) = scene.layers.current() {
                println!("layer: {}", layer.legend());
            }
            let mut picked_cell = None;
            if let Some((x, y)) = args.cursor {
                // The centre of the pixel, where the GPU samples it.
                let picked = Instant::now();
                let hit = scene.pick(
                    &view,
                    size.0,
                    size.1,
                    f64::from(x) + 0.5,
                    f64::from(y) + 0.5,
                );
                let pick_ms = picked.elapsed().as_secs_f64() * 1e3;
                match hit {
                    Some(hit) => {
                        println!(
                            "pixel {x},{y}: cell {} at {:.6} radii from the camera, picked in \
                             {pick_ms:.3} ms",
                            hit.cell, hit.distance
                        );
                        for line in inspect::cell_text(&scene.world, &scene.layers, hit.cell) {
                            println!("{line}");
                        }
                        picked_cell = Some(hit.cell);
                    }
                    None => println!("pixel {x},{y}: no cell, picked in {pick_ms:.3} ms"),
                }
            }
            let mut offscreen =
                Offscreen::new(&scene.mesh, &scene.vertex_colours(), size.0, size.1)?;
            if let Some(cell) = picked_cell {
                // Through the same write as the window's highlight.
                let colour = layers::highlight(scene.layers.cell_colour(&scene.world, cell));
                offscreen
                    .renderer()
                    .set_cell_colour(scene.mesh.cell_vertices(cell), colour);
            }
            offscreen.draw(&view.frame(aspect, &scene.relief))?;
            let mut pixels = offscreen.pixels()?;
            if let Some((x, y)) = args.cursor {
                crosshair(&mut pixels, size.0, size.1, x, y);
            }
            let png = kuniumi_io::png::encode_rgb8(size.0, size.1, &pixels);
            std::fs::write(path, png).map_err(|e| format!("{}: {e}", path.display()))?;
            println!(
                "wrote {} in {:.2} s",
                path.display(),
                start.elapsed().as_secs_f64()
            );
            Ok(())
        }
        Mode::Window => app::run(scene, view, size, None),
        Mode::Benchmark {
            frames,
            offscreen: false,
        } => app::run(scene, view, size, Some(*frames)),
        Mode::Benchmark {
            frames,
            offscreen: true,
        } => {
            let mut offscreen =
                Offscreen::new(&scene.mesh, &scene.vertex_colours(), size.0, size.1)?;
            // The first frame, not timed, sets up the targets.
            offscreen.draw(&view.frame(aspect, &scene.relief))?;
            let mut times_ms = Vec::with_capacity(*frames as usize);
            for _ in 0..*frames {
                view.camera.lon_rad += app::BENCHMARK_STEP_RAD;
                let drawn = Instant::now();
                offscreen.draw(&view.frame(aspect, &scene.relief))?;
                times_ms.push(drawn.elapsed().as_secs_f64() * 1e3);
            }
            let what = format!("{}x{} px offscreen, until the GPU finished", size.0, size.1);
            FrameStats::of(&times_ms).print(&what, scene.mesh.indices.len() / 3);
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use kuniumi_world::{PlanetParams, Settings};

    use super::*;

    fn world() -> World {
        let settings = Settings {
            cells: 1_000,
            allow_partial: true,
            ..Settings::default()
        };
        kuniumi_sim::generate(&PlanetParams::earth(), &settings, 3, &mut |_| {}).unwrap()
    }

    #[test]
    fn scene_heights_lie_above_the_sea_level() {
        let world = world();
        let sea_level = world.sea_level_m;
        assert!(sea_level < -100.0, "an Earth-like world: {sea_level}");
        let heights = world.field(SURFACE_HEIGHT).unwrap().values.clone();
        let scene = Scene::new("w".into(), world.clone(), None).unwrap();
        assert!(scene.relief.present);
        assert_eq!(
            scene.relief.radius_m.to_bits(),
            world.params.radius_m.to_bits()
        );
        assert_eq!(scene.layers.current().unwrap().name, SURFACE_HEIGHT);
        let mut top = 0.0;
        for (cell, h) in heights.iter().enumerate() {
            let above = h - sea_level;
            let centre = scene.mesh.vertices[scene.mesh.cell_vertex_offsets[cell] as usize];
            assert_eq!(centre[3].to_bits(), (above as f32).to_bits(), "cell {cell}");
            if above > top {
                top = above;
            }
        }
        assert!(top > 0.0);
        assert_eq!(scene.relief.max_height_m.to_bits(), top.to_bits());

        let plates = Scene::new("w".into(), world.clone(), Some("plate_id")).unwrap();
        assert_eq!(plates.layers.current().unwrap().name, "plate_id");
        let error = Scene::new("w".into(), world.clone(), Some("rain"))
            .err()
            .unwrap();
        assert!(
            error.starts_with("the world has no field `rain`"),
            "{error}"
        );

        // A radius that no generated world has (review R2-2).
        for radius_m in [-1.0, 0.0, f64::INFINITY] {
            let mut odd = world.clone();
            odd.params.radius_m = radius_m;
            let error = Scene::new("w".into(), odd, None).err().unwrap();
            assert_eq!(
                error,
                format!(
                    "the planet's radius is {radius_m} m; the globe needs a finite positive radius"
                )
            );
        }

        // Without heights: a sphere, no relief.
        let mut flat = world;
        flat.fields.retain(|f| f.name != SURFACE_HEIGHT);
        let scene = Scene::new("w".into(), flat, None).unwrap();
        assert!(!scene.relief.present);
        assert_eq!(scene.relief.max_height_m.to_bits(), 0.0_f64.to_bits());
        assert!(
            scene
                .mesh
                .vertices
                .iter()
                .all(|v| v[3].to_bits() == 0.0_f32.to_bits())
        );
    }

    #[test]
    fn the_scene_picks_what_the_view_shows() {
        let mut scene = Scene::new("w".into(), world(), None).unwrap();
        // The camera above the centroid of cell 42's first drawn triangle at
        // ×20: the centre of the view lies on it (a radial ray).
        let relief = 20.0;
        let sigma = f64::from((relief / scene.world.params.radius_m) as f32);
        let first = scene.mesh.cell_vertex_offsets[42] as usize;
        let drawn = |k: usize| {
            let v = scene.mesh.vertices[k];
            kuniumi_math::Vec3::new(f64::from(v[0]), f64::from(v[1]), f64::from(v[2]))
                * (1.0 + sigma * f64::from(v[3]))
        };
        let g = drawn(first) + drawn(first + 1) + drawn(first + 2);
        let (lat_rad, lon_rad) = kuniumi_sphere::coords::lat_lon(g);
        let view = View {
            camera: OrbitCamera {
                lat_rad,
                lon_rad,
                distance: 2.0,
            },
            relief,
            sea: true,
        };
        let hit = scene.pick(&view, 1280, 800, 640.0, 400.0).unwrap();
        assert_eq!(hit.cell, 42);
        let top = view.top(&scene.relief);
        assert!(hit.distance > 2.0 - top && hit.distance < 2.0, "{hit:?}");
        // The corner of the view is the background.
        let far = View {
            camera: OrbitCamera {
                distance: 3.0,
                ..view.camera
            },
            ..view
        };
        assert_eq!(scene.pick(&far, 1280, 800, 0.5, 0.5), None);
    }

    #[test]
    fn crosshair_by_hand() {
        const W: u32 = 32;
        let grey = [100, 100, 100];
        let mut pixels = vec![grey; (W * W) as usize];
        crosshair(&mut pixels, W, W, 10, 10);
        let at = |x: u32, y: u32| pixels[(y * W + x) as usize];
        let (white, black) = ([255, 255, 255], [0, 0, 0]);
        // The pixel and its neighbours untouched.
        for y in 9..=11 {
            for x in 9..=11 {
                assert_eq!(at(x, y), grey, "{x} {y}");
            }
        }
        // White arms 3 to 12 px away (to the left and up clipped at 0).
        for r in 3..=12 {
            assert_eq!(at(10 + r, 10), white);
            assert_eq!(at(10, 10 + r), white);
        }
        for r in 3..=10 {
            assert_eq!(at(10 - r, 10), white);
            assert_eq!(at(10, 10 - r), white);
        }
        // Black outline: sides and caps.
        for r in 2..=13 {
            assert_eq!(at(10 + r, 9), black);
            assert_eq!(at(10 + r, 11), black);
            assert_eq!(at(9, 10 + r), black);
            assert_eq!(at(11, 10 + r), black);
        }
        for (x, y) in [(12, 10), (23, 10), (10, 12), (10, 23), (8, 10), (10, 8)] {
            assert_eq!(at(x, y), black, "{x} {y}");
        }
        let white_count = pixels.iter().filter(|&&p| p == white).count();
        assert_eq!(white_count, 10 + 10 + 8 + 8);
        assert_eq!(at(30, 30), grey);
        assert_eq!(at(24, 10), grey);
        // In a corner, clipped without a panic.
        let mut small = vec![grey; 4 * 3];
        crosshair(&mut small, 4, 3, 3, 2);
        assert_eq!(small[(2 * 4 + 3) as usize], grey);
        assert_eq!(small[(2 * 4) as usize], white, "3 px to the left");
    }
}
