//! Search engines over an index (design §7): full scan, the pigeonhole two-copy walk, and
//! cost-model engine selection.

pub mod engine;
pub mod pigeonhole;
pub mod scan;
pub mod variants;

pub use engine::{CostModel, Engine, Plan, search};
pub use pigeonhole::pigeonhole;
pub use scan::scan;
pub use variants::Variants;
