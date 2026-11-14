//! The orbit camera: it circles the planet's centre at a latitude, longitude
//! and distance, looks at the centre and keeps north up.
//!
//! The planet is drawn with radius 1 (distances are in planet radii). The
//! matrices are computed in `f64` with `kuniumi_math::det` and rounded to
//! `f32` for the GPU at the end. Clip space is wgpu's: x to the right, y
//! up, depth from 0 at the near plane to 1 at the far plane.

use core::f64::consts::PI;

use kuniumi_math::det;
use kuniumi_math::units::RADIANS_PER_DEGREE;
use kuniumi_math::vec3::Vec3;

/// The vertical field of view (rad): 45°.
pub const FOV_Y_RAD: f64 = PI / 4.0;
/// The largest latitude of the camera (rad): 89°, so that north stays up.
pub const MAX_LAT_RAD: f64 = 89.0 * RADIANS_PER_DEGREE;
/// The largest distance from the centre, in radii of the highest drawn
/// point (planet radii for a flat globe).
pub const MAX_DISTANCE: f64 = 10.0;
/// The smallest clearance above the highest drawn point (planet radii):
/// about 64 km for Earth.
pub const MIN_CLEARANCE: f64 = 0.01;
/// The factor by which one wheel step scales the height above the surface.
pub const ZOOM_STEP: f64 = 1.1;

/// A 4 × 4 matrix as four columns (the layout of WGSL's `mat4x4<f32>`).
pub type Columns<T> = [[T; 4]; 4];

/// The camera.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrbitCamera {
    /// The latitude of the point under the camera (rad), within
    /// ±[`MAX_LAT_RAD`].
    pub lat_rad: f64,
    /// The longitude of the point under the camera (rad).
    pub lon_rad: f64,
    /// The distance from the planet's centre (planet radii), above the
    /// highest drawn point.
    pub distance: f64,
}

impl OrbitCamera {
    /// The camera's position (planet radii).
    pub fn eye(&self) -> Vec3 {
        let (sin_lat, cos_lat) = det::sin_cos(self.lat_rad);
        let (sin_lon, cos_lon) = det::sin_cos(self.lon_rad);
        Vec3::new(cos_lat * cos_lon, cos_lat * sin_lon, sin_lat) * self.distance
    }

    /// The camera's frame: forward (towards the centre), right and up, unit
    /// vectors, right-handed with +Z up.
    fn frame(&self) -> (Vec3, Vec3, Vec3) {
        let eye = self.eye();
        let forward = -eye / eye.length();
        // |lat| ≤ 89°, so forward is never parallel to +Z.
        let right = forward.cross(Vec3::Z);
        let right = right / right.length();
        (forward, right, right.cross(forward))
    }

    /// The view-projection matrix for a viewport of width/height `aspect`,
    /// for a scene whose highest point lies `top ≥ 1` radii from the centre
    /// (the near and far planes enclose the sphere of radius `top`).
    pub fn view_proj(&self, aspect: f64, top: f64) -> Columns<f64> {
        let eye = self.eye();
        let (f, s, u) = self.frame();
        // Right-handed look-at: the camera looks along −z of view space.
        let view = [
            [s.x, s.y, s.z, -s.dot(eye)],
            [u.x, u.y, u.z, -u.dot(eye)],
            [-f.x, -f.y, -f.z, f.dot(eye)],
            [0.0, 0.0, 0.0, 1.0],
        ];
        // Perspective with depth 0 at `near` and 1 at `far`; the planes
        // enclose the sphere of radius `top` (half the clearance in front).
        let near = f64::max(0.5 * (self.distance - top), 1e-4);
        let far = self.distance + top;
        let h = 1.0 / det::tan(0.5 * FOV_Y_RAD);
        let projection = [
            [h / aspect, 0.0, 0.0, 0.0],
            [0.0, h, 0.0, 0.0],
            [0.0, 0.0, far / (near - far), near * far / (near - far)],
            [0.0, 0.0, -1.0, 0.0],
        ];
        // Column `c` of projection · view, both written by rows.
        let mut columns = [[0.0; 4]; 4];
        for (c, column) in columns.iter_mut().enumerate() {
            for (r, value) in column.iter_mut().enumerate() {
                *value = (0..4).map(|k| projection[r][k] * view[k][c]).sum();
            }
        }
        columns
    }

    /// The ray from the eye through the point `ndc` of the view (normalised
    /// device coordinates: x to the right, y up, both from −1 to 1 across
    /// the viewport) of width/height `aspect`: its origin, the eye (planet
    /// radii), and its unit direction. The inverse of the perspective of
    /// [`view_proj`](Self::view_proj): `normalise(f + s·x·a/k + u·y/k)` with
    /// the forward, right and up vectors `f, s, u`, `a = aspect` and
    /// `k = 1/tan(FOV_Y/2)`.
    pub fn ray(&self, aspect: f64, ndc: [f64; 2]) -> (Vec3, Vec3) {
        let (f, s, u) = self.frame();
        let k = 1.0 / det::tan(0.5 * FOV_Y_RAD);
        let dir = f + s * (ndc[0] * aspect / k) + u * (ndc[1] / k);
        (self.eye(), dir / dir.length())
    }

    /// The direction towards the light, a unit vector: from above and to
    /// the left of the camera, so that relief facing up-left is lit.
    pub fn light_dir(&self) -> Vec3 {
        let (f, s, u) = self.frame();
        let light = -f + u * 0.6 - s * 0.4;
        light / light.length()
    }

    /// Turns the globe with a mouse drag of `(dx_px, dy_px)` (x to the
    /// right, y down) in a viewport `height_px` high: the point under the
    /// centre of the view follows the cursor (the angle per pixel is the
    /// field of view per pixel times the height above the surface).
    pub fn drag(&mut self, dx_px: f64, dy_px: f64, height_px: f64) {
        let rad_per_px = FOV_Y_RAD / height_px * (self.distance - 1.0);
        self.lon_rad -= dx_px * rad_per_px;
        self.lat_rad = (self.lat_rad + dy_px * rad_per_px).clamp(-MAX_LAT_RAD, MAX_LAT_RAD);
    }

    /// Zooms by `steps` wheel steps (positive is closer): the height above
    /// the surface is divided by [`ZOOM_STEP`] per step, then the distance
    /// is clamped ([`clamp`](Self::clamp)).
    pub fn zoom(&mut self, steps: f64, top: f64) {
        self.distance = 1.0 + (self.distance - 1.0) / det::pow(ZOOM_STEP, steps);
        self.clamp(top);
    }

    /// Clamps the distance between `top + MIN_CLEARANCE` and
    /// `MAX_DISTANCE · top` (`top ≥ 1` as for [`view_proj`](Self::view_proj),
    /// so the range is never empty) and the latitude to ±[`MAX_LAT_RAD`].
    pub fn clamp(&mut self, top: f64) {
        self.distance = self.distance.clamp(top + MIN_CLEARANCE, MAX_DISTANCE * top);
        self.lat_rad = self.lat_rad.clamp(-MAX_LAT_RAD, MAX_LAT_RAD);
    }
}

/// The normalised device coordinates `[2x/W − 1, 1 − 2y/H]` of the point
/// `(x_px, y_px)` of a viewport `width_px` × `height_px` (pixel coordinates
/// from its top-left corner, y down, as winit reports the cursor; the centre
/// of pixel `(i, j)` is `(i + ½, j + ½)`).
pub fn cursor_ndc(x_px: f64, y_px: f64, width_px: u32, height_px: u32) -> [f64; 2] {
    [
        2.0 * x_px / f64::from(width_px) - 1.0,
        1.0 - 2.0 * y_px / f64::from(height_px),
    ]
}

/// `m` rounded to `f32`.
pub fn to_f32(m: &Columns<f64>) -> Columns<f32> {
    m.map(|column| column.map(|x| x as f32))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The normalised device coordinates of `p` under `m`.
    fn project(m: &Columns<f64>, p: Vec3) -> [f64; 3] {
        let v = [p.x, p.y, p.z, 1.0];
        let clip: Vec<f64> = (0..4)
            .map(|row| (0..4).map(|col| m[col][row] * v[col]).sum())
            .collect();
        assert!(clip[3] > 0.0, "{p:?} is behind the camera");
        [clip[0] / clip[3], clip[1] / clip[3], clip[2] / clip[3]]
    }

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    const START: OrbitCamera = OrbitCamera {
        lat_rad: 0.0,
        lon_rad: 0.0,
        distance: 3.0,
    };

    #[test]
    fn values_by_hand() {
        // Camera at (3, 0, 0), looking at the centre with +Z up, so +Y (east)
        // is to the right. With top = 1 the near plane lies at half the
        // clearance, 1, and the far plane at 3 + 1 = 4.
        assert_eq!(START.eye(), Vec3::new(3.0, 0.0, 0.0));
        let m = START.view_proj(1.0, 1.0);
        let centre = project(&m, Vec3::ZERO);
        assert!(close(centre[0], 0.0) && close(centre[1], 0.0));
        // The point under the camera, 2 in front of it: depth
        // (f/(n − f)·(−2) + n·f/(n − f))/2 = (8/3 − 4/3)/2 = 2/3.
        let below = project(&m, Vec3::X);
        assert!(close(below[0], 0.0) && close(below[1], 0.0));
        assert!(close(below[2], 2.0 / 3.0), "{below:?}");
        assert!(below[2] < centre[2]);
        // East on the equator, 3 in front and 1 to the right:
        // x = cot(π/8)/3 = (1 + √2)/3.
        let east = project(&m, Vec3::Y);
        assert!(close(east[0], (1.0 + 2.0_f64.sqrt()) / 3.0), "{east:?}");
        assert!(close(east[1], 0.0));
        // North up, at the same height above the centre as east is right of it.
        let north = project(&m, Vec3::Z);
        assert!(close(north[0], 0.0) && close(north[1], east[0]));
        // The aspect squeezes x only.
        let wide = project(&START.view_proj(2.0, 1.0), Vec3::Y);
        assert!(close(wide[0], east[0] / 2.0));
        // The far side of the sphere of radius `top` lies within the depth range.
        let far = project(&m, -Vec3::X);
        assert!(far[2] > 0.0 && far[2] <= 1.0);
        let matrix = to_f32(&m);
        assert_eq!(matrix[0][0].to_bits(), (m[0][0] as f32).to_bits());
    }

    #[test]
    fn rays_project_back_to_their_screen_point() {
        let cameras = [
            START,
            OrbitCamera {
                lat_rad: 0.5,
                lon_rad: -2.0,
                distance: 4.0,
            },
            OrbitCamera {
                lat_rad: -1.5,
                lon_rad: 3.0,
                distance: 1.05,
            },
        ];
        for camera in cameras {
            for aspect in [1.0, 1.6, 0.5] {
                let m = camera.view_proj(aspect, 1.0);
                for ndc in [[0.0, 0.0], [1.0, 1.0], [-1.0, 0.25], [0.3, -0.9]] {
                    let (origin, dir) = camera.ray(aspect, ndc);
                    assert_eq!(origin, camera.eye());
                    assert!(close(dir.length(), 1.0));
                    for t in [0.5, 2.0, 5.0] {
                        let p = project(&m, origin + dir * t);
                        assert!(
                            close(p[0], ndc[0]) && close(p[1], ndc[1]),
                            "{camera:?} {aspect} {ndc:?} {t}: {p:?}"
                        );
                    }
                }
            }
        }
        // The ray through the centre of the view points at the planet's centre.
        let (origin, dir) = START.ray(1.6, [0.0, 0.0]);
        assert!((dir + origin / origin.length()).length() < 1e-15);
    }

    #[test]
    fn cursor_positions_by_hand() {
        let bits = |ndc: [f64; 2]| ndc.map(f64::to_bits);
        let at = |x, y| bits(cursor_ndc(x, y, 1280, 800));
        assert_eq!(at(0.0, 0.0), bits([-1.0, 1.0]));
        assert_eq!(at(1280.0, 800.0), bits([1.0, -1.0]));
        assert_eq!(at(640.0, 400.0), bits([0.0, 0.0]));
        assert_eq!(at(320.0, 100.0), bits([-0.5, 0.75]));
    }

    #[test]
    fn eye_and_light_follow_the_angles() {
        let camera = OrbitCamera {
            lat_rad: 0.5,
            lon_rad: -2.0,
            distance: 4.0,
        };
        let eye = camera.eye();
        let expected = Vec3::new(
            4.0 * det::cos(0.5) * det::cos(-2.0),
            4.0 * det::cos(0.5) * det::sin(-2.0),
            4.0 * det::sin(0.5),
        );
        assert!((eye - expected).length() < 1e-15);
        let light = camera.light_dir();
        assert!(close(light.length(), 1.0));
        // Mostly towards the camera, a little up (north) and to the left (west).
        assert!(light.dot(eye) > 0.0);
        let m = START.view_proj(1.0, 1.0);
        let lit = project(&m, START.light_dir());
        assert!(lit[0] < 0.0 && lit[1] > 0.0, "{lit:?}");
    }

    #[test]
    fn dragging_moves_the_globe_with_the_cursor() {
        let height_px = 800.0;
        let mut right = START;
        right.drag(100.0, 0.0, height_px);
        // Height above the surface 2: 100 px of a 45° view over 800 px turn
        // the globe by 100 · (π/4)/800 · 2 rad.
        let angle = 100.0 * FOV_Y_RAD / height_px * 2.0;
        assert!(close(right.lon_rad, -angle) && right.lat_rad == 0.0);
        let p = project(&right.view_proj(1.0, 1.0), Vec3::X);
        assert!(
            p[0] > 0.0 && close(p[1], 0.0),
            "the point moved right: {p:?}"
        );
        let mut down = START;
        down.drag(0.0, 100.0, height_px);
        assert!(close(down.lat_rad, angle) && down.lon_rad == 0.0);
        let p = project(&down.view_proj(1.0, 1.0), Vec3::X);
        assert!(
            p[1] < 0.0 && close(p[0], 0.0),
            "the point moved down: {p:?}"
        );
        // Closer to the surface, slower.
        let mut near = OrbitCamera {
            distance: 1.5,
            ..START
        };
        near.drag(100.0, 0.0, height_px);
        assert!(close(near.lon_rad, -angle / 4.0));
        // The latitude stops at ±89°.
        let mut pole = START;
        pole.drag(0.0, 1e6, height_px);
        assert_eq!(pole.lat_rad.to_bits(), MAX_LAT_RAD.to_bits());
        pole.drag(0.0, -1e7, height_px);
        assert_eq!(pole.lat_rad.to_bits(), (-MAX_LAT_RAD).to_bits());
    }

    #[test]
    fn zoom_scales_the_height_and_clamps() {
        let mut camera = START;
        camera.zoom(1.0, 1.0);
        assert!(close(camera.distance, 1.0 + 2.0 / 1.1));
        camera.zoom(-2.0, 1.0);
        assert!(
            close(camera.distance, 1.0 + 2.0 * 1.1),
            "{}",
            camera.distance
        );
        camera.zoom(1_000.0, 1.2);
        assert_eq!(camera.distance.to_bits(), (1.2 + MIN_CLEARANCE).to_bits());
        camera.zoom(-1_000.0, 1.2);
        assert_eq!(camera.distance.to_bits(), (MAX_DISTANCE * 1.2).to_bits());
        camera.zoom(-1_000.0, 1.0);
        assert_eq!(camera.distance.to_bits(), MAX_DISTANCE.to_bits());
        // A higher relief pushes the camera out.
        let mut low = OrbitCamera {
            distance: 1.05,
            ..START
        };
        low.clamp(1.1);
        assert_eq!(low.distance.to_bits(), (1.1 + MIN_CLEARANCE).to_bits());
        let mut tilted = OrbitCamera {
            lat_rad: 2.0,
            ..START
        };
        tilted.clamp(1.0);
        assert_eq!(tilted.lat_rad.to_bits(), MAX_LAT_RAD.to_bits());
        // Relief far above 10 radii (review R1-1): the range grows with it.
        let mut tall = START;
        tall.clamp(20.0);
        assert_eq!(tall.distance.to_bits(), (20.0 + MIN_CLEARANCE).to_bits());
        tall.zoom(-1_000.0, 20.0);
        assert_eq!(tall.distance.to_bits(), (MAX_DISTANCE * 20.0).to_bits());
    }
}
