//! File formats and image output.
//!
//! - [`png`]: a dependency-free PNG encoder (uncompressed deflate), used for
//!   debug images that both humans and Claude Code can look at.
//! - [`field_file`]: the `.kfield` binary format for one per-cell field, and a
//!   bit-exact fingerprint used by reproducibility tests.
//!
//! Everything that touches the file system lives in this crate (and the CLI),
//! never in `kuniumi-sim`.

pub mod field_file;
pub mod png;
