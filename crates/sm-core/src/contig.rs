//! Global coordinate space over all FASTA records (design §5.2).

/// One FASTA record placed in the global coordinate space.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Contig {
    /// Record name: the first word of the FASTA header.
    pub name: String,
    /// Rest of the header line, if any.
    pub description: String,
    /// Global coordinate of the record's first base.
    pub offset: u32,
    pub len: u32,
}

impl Contig {
    pub fn end(&self) -> u32 {
        self.offset + self.len
    }
}

/// Contig table: records in FASTA order, laid end to end with no separators.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Contigs {
    contigs: Vec<Contig>,
}

#[derive(Debug, thiserror::Error)]
#[error("genome exceeds {max} bases (the u32 coordinate limit)", max = u32::MAX)]
pub struct TooLarge;

impl Contigs {
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a record of `len` bases; fails if the genome would exceed `u32::MAX`.
    pub fn push(&mut self, name: String, description: String, len: u64) -> Result<&Contig, TooLarge> {
        let offset = self.total_len();
        let end = offset + len;
        if end > u32::MAX as u64 {
            return Err(TooLarge);
        }
        self.contigs.push(Contig { name, description, offset: offset as u32, len: len as u32 });
        Ok(self.contigs.last().unwrap())
    }

    pub fn total_len(&self) -> u64 {
        self.contigs.last().map_or(0, |c| c.end() as u64)
    }

    pub fn len(&self) -> usize {
        self.contigs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.contigs.is_empty()
    }

    pub fn iter(&self) -> std::slice::Iter<'_, Contig> {
        self.contigs.iter()
    }

    pub fn get(&self, i: usize) -> &Contig {
        &self.contigs[i]
    }

    pub fn by_name(&self, name: &str) -> Option<usize> {
        self.contigs.iter().position(|c| c.name == name)
    }

    /// Index of the contig containing global position `pos`.
    /// Zero-length contigs are never returned.
    pub fn locate(&self, pos: u32) -> Option<usize> {
        let i = self.contigs.partition_point(|c| c.end() <= pos);
        (i < self.contigs.len() && self.contigs[i].offset <= pos).then_some(i)
    }

    /// `(contig index, 0-based offset within it)` for a global position.
    pub fn to_local(&self, pos: u32) -> Option<(usize, u32)> {
        self.locate(pos).map(|i| (i, pos - self.contigs[i].offset))
    }
}

impl<'a> IntoIterator for &'a Contigs {
    type Item = &'a Contig;
    type IntoIter = std::slice::Iter<'a, Contig>;
    fn into_iter(self) -> Self::IntoIter {
        self.contigs.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn locate() {
        let mut c = Contigs::new();
        c.push("a".into(), String::new(), 10).unwrap();
        c.push("empty".into(), String::new(), 0).unwrap();
        c.push("b".into(), String::new(), 5).unwrap();
        assert_eq!(c.to_local(0), Some((0, 0)));
        assert_eq!(c.to_local(9), Some((0, 9)));
        assert_eq!(c.to_local(10), Some((2, 0)));
        assert_eq!(c.to_local(14), Some((2, 4)));
        assert_eq!(c.to_local(15), None);
        assert_eq!(c.total_len(), 15);
    }

    #[test]
    fn rejects_oversize() {
        let mut c = Contigs::new();
        c.push("a".into(), String::new(), u32::MAX as u64).unwrap();
        assert!(c.push("b".into(), String::new(), 1).is_err());
    }
}
