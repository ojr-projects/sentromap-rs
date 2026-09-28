//! On-disk index formats (design §6), memory-mapped readers and the partitioned builder (§13).

pub mod array;
pub mod bits;
pub mod build;
pub mod index;
pub mod kmers;
pub mod manifest;
pub mod positions;
pub mod seq;
pub mod verify;

pub use build::{BuildOptions, build};
pub use index::Index;
pub use manifest::Manifest;
