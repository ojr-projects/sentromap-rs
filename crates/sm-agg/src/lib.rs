//! Server-side result handling: result sets (design §8), orderings per n (§9), and aggregation
//! for display (§10): screen-resolution rasters and per-bin histograms.

pub mod order;
pub mod raster;
pub mod result;

pub use order::{OrderMode, OrderView};
pub use raster::{Raster, View, histogram, raster};
pub use result::{ResultSet, SiteTable};
