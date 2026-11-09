//! File formats and image output.
//!
//! - [`colormap`]: colour maps from values to 8-bit sRGB ([`Rgb8`]):
//!   viridis, a diverging map, an elevation map with a break at sea level and
//!   a categorical palette.
//! - [`render`]: equirectangular and orthographic (globe) maps of per-cell
//!   fields on a [`kuniumi_sphere::SphereGrid`].
//! - [`png`]: compressed PNG encoding (8-bit RGB and 16-bit grey), through
//!   the `png` crate, for images that both humans and Claude Code look at.
//! - [`field_file`]: the `.kfield` binary format for one per-cell field, and a
//!   bit-exact fingerprint used by reproducibility tests.
//! - [`world_dir`]: a world as a directory: a `.kfield` file and a PNG map
//!   per field, and `world.json` with the metadata (through `serde_json`).
//!
//! Everything that touches the file system lives in this crate (and the CLI),
//! never in `kuniumi-sim`.

pub mod colormap;
pub mod field_file;
pub mod png;
pub mod render;
pub mod world_dir;

pub use colormap::Rgb8;
