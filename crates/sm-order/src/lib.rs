//! Ordering search variants as a cladogram (design §9): the shared-character tree for any
//! size, and an exact MST for small sets.

pub mod chartree;
pub mod mst;

pub use chartree::{Clade, Ordered, Ranking, filter_frozen, order, order_per_n, star_length};
pub use mst::{MST_LIMIT, Mst, mst};
