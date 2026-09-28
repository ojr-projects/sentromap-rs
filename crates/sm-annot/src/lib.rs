//! Genome annotations (design §11): interval features (GFF3, RepeatMasker, BED) and
//! quantitative tracks (bigWig, GC%), imported into the index directory once and memory-mapped
//! by the server.

pub mod alias;
pub mod features;
pub mod tracks;

use std::path::Path;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

pub use alias::Aliases;
pub use features::{FeatureSet, FeatureSetInfo};
pub use tracks::{Track, TrackInfo};

pub const REGISTRY: &str = "annotations.json";

/// What has been imported into an index directory.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Registry {
    /// Short display name per contig (from alias tables).
    pub contig_display: Vec<String>,
    pub feature_sets: Vec<FeatureSetInfo>,
    pub tracks: Vec<TrackInfo>,
}

impl Registry {
    pub fn load(dir: &Path) -> Result<Self> {
        let p = dir.join(REGISTRY);
        if !p.exists() {
            return Ok(Self::default());
        }
        let s = std::fs::read_to_string(&p)?;
        serde_json::from_str(&s).with_context(|| format!("parsing {}", p.display()))
    }

    pub fn save(&self, dir: &Path) -> Result<()> {
        sm_index::manifest::write_json(&dir.join(REGISTRY), self)
    }

    pub fn upsert_features(&mut self, info: FeatureSetInfo) {
        self.feature_sets.retain(|f| f.name != info.name);
        self.feature_sets.push(info);
    }

    pub fn upsert_track(&mut self, info: TrackInfo) {
        self.tracks.retain(|t| t.name != info.name);
        self.tracks.push(info);
    }
}

/// Opened annotations.
pub struct Annotations {
    pub registry: Registry,
    pub feature_sets: Vec<FeatureSet>,
    pub tracks: Vec<Track>,
}

impl Annotations {
    pub fn open(dir: &Path) -> Result<Self> {
        let registry = Registry::load(dir)?;
        let feature_sets =
            registry.feature_sets.iter().map(|i| FeatureSet::open(dir, i.clone())).collect::<Result<_>>()?;
        let tracks = registry.tracks.iter().map(|i| Track::open(dir, i.clone())).collect::<Result<_>>()?;
        Ok(Self { registry, feature_sets, tracks })
    }

    pub fn feature_set(&self, name: &str) -> Option<&FeatureSet> {
        self.feature_sets.iter().find(|f| f.info.name == name)
    }

    pub fn track(&self, name: &str) -> Option<&Track> {
        self.tracks.iter().find(|t| t.info.name == name)
    }
}
