//! `kuniumi-viewer`: a world directory as an interactive 3D globe (ticket
//! MV-02; D-19).
//!
//! `kuniumi-viewer <world dir>` reads the world (`kuniumi_io::world_dir`),
//! builds a triangle mesh of its Voronoi cells ([`mesh`]) and shows any of
//! its fields ([`layers`]) in the colours of the world's PNG maps, with the
//! relief exaggerated and a translucent sea-level surface ([`view`]), on the
//! GPU through `wgpu` ([`gpu`]) in a `winit` window ([`app`]). The mouse
//! turns and zooms the globe ([`camera`]), the keyboard changes the layer,
//! the relief and the sea. `--screenshot` renders one frame into a PNG
//! without a window; `--benchmark` times the frames of a spinning globe.
//! The usage is in [`args::USAGE`].

// The GPU works in `f32`. D-19 and rule 1 of docs/conventions.md,
// "Determinism", let this crate (and kuniumi-gpu) use `f32` through this
// crate-level allowance: nothing the viewer computes feeds back into
// generation, and no other crate may do this.
#![allow(clippy::disallowed_types)]

mod app;
mod args;
mod camera;
mod gpu;
mod layers;
mod mesh;
mod view;

use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Instant;

use kuniumi_sim::stages::crust::SURFACE_HEIGHT;
use kuniumi_world::World;

use crate::app::FrameStats;
use crate::args::{Args, Mode, Parsed, USAGE};
use crate::camera::OrbitCamera;
use crate::gpu::Offscreen;
use crate::layers::Layers;
use crate::mesh::GlobeMesh;
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
        Ok(Self {
            dir,
            world,
            mesh,
            layers,
            relief,
        })
    }

    /// The colour of every vertex in the layer shown.
    pub fn vertex_colours(&self) -> Vec<[u8; 4]> {
        self.mesh.vertex_colours(&self.layers.colours(&self.world))
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
    let scene = Scene::new(args.world_dir.clone(), world, args.layer.as_deref())?;
    println!(
        "{}: {} cells, {} fields; read in {read_s:.2} s, mesh of {} triangles and layers in {:.2} s",
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
            let mut offscreen =
                Offscreen::new(&scene.mesh, &scene.vertex_colours(), size.0, size.1)?;
            offscreen.draw(&view.frame(aspect, &scene.relief))?;
            let png = kuniumi_io::png::encode_rgb8(size.0, size.1, &offscreen.pixels()?);
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
}
