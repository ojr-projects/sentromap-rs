//! Chromosome name aliases: map UCSC/Ensembl/GenBank/RefSeq names onto the index's contigs.

use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};
use sm_core::Contigs;

pub struct Aliases {
    to_contig: HashMap<String, usize>,
    /// A short display name per contig (the shortest alias, e.g. "2L" for NT_033779.5).
    pub display: Vec<String>,
}

impl Aliases {
    /// Contig names only.
    pub fn new(contigs: &Contigs) -> Self {
        let to_contig = contigs.iter().enumerate().map(|(i, c)| (c.name.clone(), i)).collect();
        let display = contigs.iter().map(|c| c.name.clone()).collect();
        Self { to_contig, display }
    }

    /// Add a tab-separated alias table (UCSC `chromAlias.txt` style): every row lists names of
    /// one sequence; a row naming a contig maps all its names to that contig.
    pub fn load(&mut self, path: &Path) -> Result<usize> {
        let text = std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
        let mut rows = 0;
        for line in text.lines() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let names: Vec<&str> = line.split('\t').map(str::trim).filter(|s| !s.is_empty()).collect();
            let Some(&ci) = names.iter().find_map(|n| self.to_contig.get(*n)) else { continue };
            for n in &names {
                self.to_contig.entry(n.to_string()).or_insert(ci);
            }
            if let Some(short) = names.iter().min_by_key(|n| (n.len(), n.to_string()))
                && short.len() < self.display[ci].len()
            {
                self.display[ci] = short.to_string();
            }
            rows += 1;
        }
        Ok(rows)
    }

    pub fn contig(&self, name: &str) -> Option<usize> {
        self.to_contig.get(name).copied()
    }
}
