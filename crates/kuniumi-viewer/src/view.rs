//! What the viewer shows besides the layer: the camera, the relief
//! exaggeration and the sea-level surface, and the values a frame hands to
//! the GPU.
//!
//! A vertex at height `h` above the sea level (m) is drawn at the radius
//! `1 + e·h/R` (planet radii), with the exaggeration `e` and the planet's
//! radius `R` (m); the sea surface is the same mesh at `h = 0`, radius 1.
//! Each world caps `e` so that its deepest vertex is drawn no closer to the
//! centre than [`MIN_DRAWN_RADIUS`] ([`Relief::max_relief`], ticket MV-05).

use crate::camera::{Columns, OrbitCamera, to_f32};

/// The relief exaggerations that Up and Down step through.
pub const RELIEF_STEPS: [f64; 8] = [0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0];

/// The largest relief exaggeration on any world, the largest key step: at
/// ×1000, Earth's deepest point, about 11 km below the sea level (the
/// Challenger Deep, 10,935 m) on a radius of 6,371 km, would be drawn
/// beyond the centre (1 − 1000 × 10.9/6,371 = −0.72), inside out. Each
/// world lowers it further ([`Relief::max_relief`]).
pub const MAX_RELIEF: f64 = RELIEF_STEPS[RELIEF_STEPS.len() - 1];

/// The smallest radius at which a vertex is drawn (planet radii): the
/// relief exaggeration of each world is capped so that its deepest vertex
/// stays at half the radius or above ([`Relief::max_relief`]).
pub const MIN_DRAWN_RADIUS: f64 = 0.5;

/// The relative margin by which [`Relief::max_relief`] stays below the
/// exaggeration that would draw the deepest vertex exactly at
/// [`MIN_DRAWN_RADIUS`]: room for the GPU's `f32` arithmetic (see there).
pub const RELIEF_MARGIN: f64 = 1e-6;

/// The smallest step above `relief`, lowered to `max`; `relief` itself if no
/// step lies above it or it is already at `max` or above.
pub fn relief_up(relief: f64, max: f64) -> f64 {
    match RELIEF_STEPS.into_iter().find(|&step| step > relief) {
        Some(step) if relief < max => step.min(max),
        _ => relief,
    }
}

/// The largest step below `relief`, or 0.
pub fn relief_down(relief: f64) -> f64 {
    RELIEF_STEPS
        .into_iter()
        .rev()
        .find(|&step| step < relief)
        .unwrap_or(0.0)
}

/// `relief` as the title and the terminal show it: rounded to two decimals,
/// without trailing zeros (`20`, `14.8`).
pub fn relief_text(relief: f64) -> String {
    format!("{}", (relief * 100.0).round() / 100.0)
}

/// The terminal line after Up or Down set the relief to `relief` on a world
/// whose largest is `max` ([`Relief::max_relief`]): the relief, and at a
/// cap below [`MAX_RELIEF`] that it is the most for this world.
pub fn relief_line(relief: f64, max: f64) -> String {
    if relief >= max && max < MAX_RELIEF {
        format!(
            "relief x{}, the most for this world (its deepest point at {} of the radius)",
            relief_text(relief),
            MIN_DRAWN_RADIUS
        )
    } else {
        format!("relief x{}", relief_text(relief))
    }
}

/// The relief of the world as drawn: the planet radius and the heights.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Relief {
    /// The planet's radius (m).
    pub radius_m: f64,
    /// The largest height above the sea level (m), at least 0.
    pub max_height_m: f64,
    /// The lowest height above the sea level (m), at most 0: minus the depth
    /// of the deepest point below the sea level.
    pub min_height_m: f64,
    /// Whether the world has heights (`surface_height`); without them the
    /// globe is a sphere and there is no sea surface.
    pub present: bool,
}

impl Relief {
    /// The largest relief exaggeration for this world: [`MAX_RELIEF`], or
    /// less where that would draw the deepest point closer to the centre
    /// than [`MIN_DRAWN_RADIUS`].
    ///
    /// With the depth `d = −min_height_m > 0` of the deepest point below the
    /// sea level and the radius `R`, `e_max = min(MAX_RELIEF,
    /// (1 − MIN_DRAWN_RADIUS)·R/d·(1 − δ))`, `δ` = [`RELIEF_MARGIN`]; a world
    /// without heights or without a point below the sea level gets
    /// [`MAX_RELIEF`]. At `e ≤ e_max` every vertex is drawn at a radius
    /// `1 + e·h/R ≥ 1 − e·d/R ≥ MIN_DRAWN_RADIUS + δ/2` (h ≥ −d). The GPU
    /// computes the radius in `f32` from `σ = fl32(e/R)` and the `f32`
    /// height: the roundings of `σ`, of `h`, of the product and of the sum
    /// (at most 2⁻²⁴ relative each, or one rounding with a fused
    /// multiply-add) move it near 0.5 by less than 1.2·10⁻⁷, below
    /// `δ/2` = 5·10⁻⁷, so the drawn radius is at least
    /// [`MIN_DRAWN_RADIUS`].
    ///
    /// By hand: a planet of 1,500 km whose sea floor reaches 50.7 km below
    /// the sea level gets 0.5 × 1,500/50.7 = 14.8; Earth, 10.9 km deep at
    /// the Challenger Deep, 0.5 × 6,371/10.9 = 291, so [`MAX_RELIEF`].
    pub fn max_relief(&self) -> f64 {
        let depth_m = -self.min_height_m;
        if !self.present || depth_m <= 0.0 {
            return MAX_RELIEF;
        }
        let cap = (1.0 - MIN_DRAWN_RADIUS) * self.radius_m / depth_m * (1.0 - RELIEF_MARGIN);
        cap.min(MAX_RELIEF)
    }
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
        let up = |relief: f64| relief_up(relief, MAX_RELIEF);
        assert_eq!(up(0.0).to_bits(), 1.0_f64.to_bits());
        assert_eq!(up(20.0).to_bits(), 50.0_f64.to_bits());
        assert_eq!(up(3.0).to_bits(), 5.0_f64.to_bits());
        assert_eq!(up(100.0).to_bits(), 100.0_f64.to_bits());
        assert_eq!(up(150.0).to_bits(), 150.0_f64.to_bits());
        assert_eq!(relief_down(20.0).to_bits(), 10.0_f64.to_bits());
        assert_eq!(relief_down(3.0).to_bits(), 2.0_f64.to_bits());
        assert_eq!(relief_down(0.5).to_bits(), 0.0_f64.to_bits());
        assert_eq!(relief_down(0.0).to_bits(), 0.0_f64.to_bits());
        assert_eq!(relief_down(150.0).to_bits(), 100.0_f64.to_bits());
        // Under a world's cap (MV-05): the next step lowered to the cap, and
        // no step beyond it; Up never lowers the relief.
        assert_eq!(relief_up(10.0, 14.8).to_bits(), 14.8_f64.to_bits());
        assert_eq!(relief_up(14.8, 14.8).to_bits(), 14.8_f64.to_bits());
        assert_eq!(relief_up(5.0, 14.8).to_bits(), 10.0_f64.to_bits());
        assert_eq!(relief_up(0.0, 0.5).to_bits(), 0.5_f64.to_bits());
        assert_eq!(relief_up(0.5, 0.5).to_bits(), 0.5_f64.to_bits());
        assert_eq!(relief_up(20.0, 14.8).to_bits(), 20.0_f64.to_bits());
        assert_eq!(relief_down(14.8).to_bits(), 10.0_f64.to_bits());
        assert_eq!(relief_down(0.25).to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn relief_limits_by_hand() {
        assert_eq!(MAX_RELIEF.to_bits(), 100.0_f64.to_bits());
        // 1,500 km with the deepest point 50 km below the sea level:
        // (1 − 0.5) × 1,500/50 = 15, less the margin; drawn at 0.5 + 5e-7.
        let small = Relief {
            radius_m: 1_500_000.0,
            max_height_m: 0.0,
            min_height_m: -50_000.0,
            present: true,
        };
        let cap = small.max_relief();
        assert!((cap - 15.0 * (1.0 - 1e-6)).abs() < 1e-12, "{cap}");
        assert!(cap < 15.0);
        let radius = 1.0 + cap / small.radius_m * small.min_height_m;
        assert!((radius - (0.5 + 5e-7)).abs() < 1e-12, "{radius}");
        // Earth with the Challenger Deep, 10,935 m below the sea level:
        // 0.5 × 6,371/10.935 = 291, above the largest step.
        let earth = Relief {
            radius_m: 6_371_000.0,
            max_height_m: 8_849.0,
            min_height_m: -10_935.0,
            present: true,
        };
        assert_eq!(earth.max_relief().to_bits(), 100.0_f64.to_bits());
        // Just deeper than what ×100 allows on Earth's radius:
        // 0.5 × 6,371/32 = 99.55.
        let deeper = Relief {
            min_height_m: -32_000.0,
            ..earth
        };
        let cap = deeper.max_relief();
        assert!((cap - 99.546_875 * (1.0 - 1e-6)).abs() < 1e-11, "{cap}");
        // A floor deeper than the radius: below ×1, 0.5 × 1,500/3,000.
        let abyss = Relief {
            min_height_m: -3_000_000.0,
            ..small
        };
        let cap = abyss.max_relief();
        assert!((cap - 0.25 * (1.0 - 1e-6)).abs() < 1e-15, "{cap}");
        // Nothing below the sea level, or no heights: the largest step.
        let land = Relief {
            min_height_m: 0.0,
            ..small
        };
        assert_eq!(land.max_relief().to_bits(), 100.0_f64.to_bits());
        let sphere = Relief {
            present: false,
            ..small
        };
        assert_eq!(sphere.max_relief().to_bits(), 100.0_f64.to_bits());
    }

    #[test]
    fn relief_lines() {
        // Review R1-1: at a world's cap, and below it or at ×100.
        let cap = 14.798_266;
        assert_eq!(
            relief_line(cap, cap),
            "relief x14.8, the most for this world (its deepest point at 0.5 of the radius)"
        );
        assert_eq!(relief_line(10.0, cap), "relief x10");
        assert_eq!(relief_line(100.0, MAX_RELIEF), "relief x100");
        assert_eq!(relief_line(20.0, MAX_RELIEF), "relief x20");
    }

    #[test]
    fn relief_texts() {
        assert_eq!(relief_text(20.0), "20");
        assert_eq!(relief_text(0.0), "0");
        assert_eq!(relief_text(100.0), "100");
        assert_eq!(relief_text(14.798_280_787), "14.8");
        assert_eq!(relief_text(0.25), "0.25");
        assert_eq!(relief_text(1.0 / 3.0), "0.33");
        assert_eq!(relief_text(99.999_9), "100");
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
            min_height_m: -11_000.0,
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
