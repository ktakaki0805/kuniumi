//! What the viewer shows besides the layer: the camera, the relief
//! exaggeration and the sea-level surface, and the values a frame hands to
//! the GPU.
//!
//! A vertex at height `h` above the sea level (m) is drawn at the radius
//! `1 + e·h/R` (planet radii), with the exaggeration `e` and the planet's
//! radius `R` (m); the sea surface is the same mesh at `h = 0`, radius 1.

use crate::camera::{Columns, OrbitCamera, to_f32};

/// The relief exaggerations that Up and Down step through.
pub const RELIEF_STEPS: [f64; 8] = [0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0];

/// The smallest step above `relief`, or `relief` itself above the last step.
pub fn relief_up(relief: f64) -> f64 {
    RELIEF_STEPS
        .into_iter()
        .find(|&step| step > relief)
        .unwrap_or(relief)
}

/// The largest step below `relief`, or 0.
pub fn relief_down(relief: f64) -> f64 {
    RELIEF_STEPS
        .into_iter()
        .rev()
        .find(|&step| step < relief)
        .unwrap_or(0.0)
}

/// The relief of the world as drawn: the planet radius and the heights.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Relief {
    /// The planet's radius (m).
    pub radius_m: f64,
    /// The largest height above the sea level (m), at least 0.
    pub max_height_m: f64,
    /// Whether the world has heights (`surface_height`); without them the
    /// globe is a sphere and there is no sea surface.
    pub present: bool,
}

/// The view: camera, relief exaggeration, sea surface on or off.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct View {
    /// The camera.
    pub camera: OrbitCamera,
    /// The relief exaggeration `e`, at least 0.
    pub relief: f64,
    /// Whether the sea-level surface is wanted.
    pub sea: bool,
}

/// The values of one frame for the GPU.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Frame {
    /// The view-projection matrix.
    pub view_proj: Columns<f32>,
    /// The direction towards the light (unit vector).
    pub light_dir: [f32; 3],
    /// `e/R`, the radius added per metre of height (planet radii per m).
    pub relief_scale_per_m: f32,
    /// Whether the sea pass runs: wanted, and the world has relief drawn
    /// with `e > 0` (at `e = 0` the sea surface would coincide with the
    /// ground).
    pub sea: bool,
}

impl View {
    /// `e/R` for `relief` (0 for a world without heights).
    pub fn relief_scale_per_m(&self, relief: &Relief) -> f64 {
        if relief.present {
            self.relief / relief.radius_m
        } else {
            0.0
        }
    }

    /// The radius of the highest drawn point (planet radii), at least 1.
    pub fn top(&self, relief: &Relief) -> f64 {
        1.0 + self.relief_scale_per_m(relief) * relief.max_height_m
    }

    /// The frame for a viewport of width/height `aspect`.
    pub fn frame(&self, aspect: f64, relief: &Relief) -> Frame {
        let scale = self.relief_scale_per_m(relief);
        let light = self.camera.light_dir();
        Frame {
            view_proj: to_f32(&self.camera.view_proj(aspect, self.top(relief))),
            light_dir: [light.x as f32, light.y as f32, light.z as f32],
            relief_scale_per_m: scale as f32,
            sea: self.sea && scale > 0.0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn relief_steps() {
        assert_eq!(relief_up(0.0).to_bits(), 1.0_f64.to_bits());
        assert_eq!(relief_up(20.0).to_bits(), 50.0_f64.to_bits());
        assert_eq!(relief_up(3.0).to_bits(), 5.0_f64.to_bits());
        assert_eq!(relief_up(100.0).to_bits(), 100.0_f64.to_bits());
        assert_eq!(relief_up(150.0).to_bits(), 150.0_f64.to_bits());
        assert_eq!(relief_down(20.0).to_bits(), 10.0_f64.to_bits());
        assert_eq!(relief_down(3.0).to_bits(), 2.0_f64.to_bits());
        assert_eq!(relief_down(0.5).to_bits(), 0.0_f64.to_bits());
        assert_eq!(relief_down(0.0).to_bits(), 0.0_f64.to_bits());
        assert_eq!(relief_down(150.0).to_bits(), 100.0_f64.to_bits());
    }

    #[test]
    fn frames() {
        let camera = OrbitCamera {
            lat_rad: 0.3,
            lon_rad: 1.0,
            distance: 2.5,
        };
        let earth = Relief {
            radius_m: 6_371_000.0,
            max_height_m: 6_371.0,
            present: true,
        };
        let view = View {
            camera,
            relief: 20.0,
            sea: true,
        };
        assert_eq!(
            view.relief_scale_per_m(&earth).to_bits(),
            (20.0_f64 / 6_371_000.0).to_bits()
        );
        // 20 × 6,371 m on a radius of 6,371 km: 2 % above the sphere.
        assert!((view.top(&earth) - 1.02).abs() < 1e-15);
        let frame = view.frame(1.5, &earth);
        assert_eq!(
            frame.view_proj,
            to_f32(&camera.view_proj(1.5, view.top(&earth)))
        );
        let light = camera.light_dir();
        assert_eq!(
            frame.light_dir.map(f32::to_bits),
            [light.x as f32, light.y as f32, light.z as f32].map(f32::to_bits)
        );
        assert_eq!(
            frame.relief_scale_per_m.to_bits(),
            ((20.0 / 6_371_000.0) as f32).to_bits()
        );
        assert!(frame.sea);
        // No sea surface when it is switched off, at e = 0 or without heights.
        assert!(!View { sea: false, ..view }.frame(1.5, &earth).sea);
        let flat = View {
            relief: 0.0,
            ..view
        };
        assert!(!flat.frame(1.5, &earth).sea);
        assert_eq!(flat.top(&earth).to_bits(), 1.0_f64.to_bits());
        let sphere = Relief {
            present: false,
            ..earth
        };
        assert_eq!(
            view.relief_scale_per_m(&sphere).to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(view.top(&sphere).to_bits(), 1.0_f64.to_bits());
        let frame = view.frame(1.5, &sphere);
        assert!(!frame.sea);
        assert_eq!(frame.relief_scale_per_m.to_bits(), 0.0_f32.to_bits());
    }
}
