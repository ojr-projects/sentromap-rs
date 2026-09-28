//! Search engines over an index: full scan (§7.3), with more to come (pigeonhole walk §7.4,
//! cost-model selection §7.5).

pub mod scan;
pub mod variants;

pub use scan::scan;
pub use variants::Variants;
