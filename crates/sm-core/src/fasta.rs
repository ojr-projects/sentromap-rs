//! Minimal streaming FASTA reader (plain or gzip) and an in-memory genome.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

use crate::contig::Contigs;

/// Open a file for buffered reading, transparently decompressing gzip.
pub fn open(path: &Path) -> io::Result<Box<dyn BufRead + Send>> {
    let mut f = File::open(path)?;
    let mut magic = [0u8; 2];
    let n = f.read(&mut magic)?;
    drop(f);
    let f = File::open(path)?;
    if n == 2 && magic == [0x1f, 0x8b] {
        Ok(Box::new(BufReader::with_capacity(1 << 20, flate2::read::MultiGzDecoder::new(f))))
    } else {
        Ok(Box::new(BufReader::with_capacity(1 << 20, f)))
    }
}

/// Header of one FASTA record.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Header {
    pub name: String,
    pub description: String,
}

/// Reads one record at a time into a caller-provided buffer.
pub struct FastaReader<R> {
    inner: R,
    line: Vec<u8>,
    pending: Option<Header>,
    started: bool,
}

impl<R: BufRead> FastaReader<R> {
    pub fn new(inner: R) -> Self {
        Self { inner, line: Vec::new(), pending: None, started: false }
    }

    /// Read the next record's sequence (raw bytes, line breaks removed) into `seq`.
    pub fn next_record(&mut self, seq: &mut Vec<u8>) -> io::Result<Option<Header>> {
        seq.clear();
        if !self.started {
            self.started = true;
            loop {
                self.line.clear();
                if self.inner.read_until(b'\n', &mut self.line)? == 0 {
                    return Ok(None);
                }
                if self.line.starts_with(b">") {
                    self.pending = Some(parse_header(&self.line));
                    break;
                }
                if !trim(&self.line).is_empty() {
                    return Err(io::Error::new(io::ErrorKind::InvalidData, "FASTA must start with '>'"));
                }
            }
        }
        let Some(header) = self.pending.take() else { return Ok(None) };
        loop {
            self.line.clear();
            if self.inner.read_until(b'\n', &mut self.line)? == 0 {
                break;
            }
            if self.line.starts_with(b">") {
                self.pending = Some(parse_header(&self.line));
                break;
            }
            seq.extend_from_slice(trim(&self.line));
        }
        Ok(Some(header))
    }
}

fn trim(line: &[u8]) -> &[u8] {
    let mut end = line.len();
    while end > 0 && line[end - 1].is_ascii_whitespace() {
        end -= 1;
    }
    &line[..end]
}

fn parse_header(line: &[u8]) -> Header {
    let text = String::from_utf8_lossy(trim(&line[1..]));
    let text = text.trim();
    match text.split_once(char::is_whitespace) {
        Some((n, d)) => Header { name: n.to_string(), description: d.trim().to_string() },
        None => Header { name: text.to_string(), description: String::new() },
    }
}

/// A whole genome in memory: every record's bases concatenated in global coordinates.
#[derive(Clone, Debug, Default)]
pub struct Genome {
    pub contigs: Contigs,
    pub seq: Vec<u8>,
}

impl Genome {
    pub fn read(path: &Path) -> io::Result<Self> {
        let mut r = FastaReader::new(open(path)?);
        let mut g = Genome::default();
        let mut buf = Vec::new();
        while let Some(h) = r.next_record(&mut buf)? {
            g.push(h.name, h.description, &buf)?;
        }
        Ok(g)
    }

    pub fn push(&mut self, name: String, description: String, seq: &[u8]) -> io::Result<()> {
        self.contigs
            .push(name, description, seq.len() as u64)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        self.seq.extend_from_slice(seq);
        Ok(())
    }

    /// The bases of contig `i`.
    pub fn contig_seq(&self, i: usize) -> &[u8] {
        let c = self.contigs.get(i);
        &self.seq[c.offset as usize..c.end() as usize]
    }

    pub fn write_fasta(&self, w: &mut impl Write) -> io::Result<()> {
        for (i, c) in self.contigs.iter().enumerate() {
            if c.description.is_empty() {
                writeln!(w, ">{}", c.name)?;
            } else {
                writeln!(w, ">{} {}", c.name, c.description)?;
            }
            for line in self.contig_seq(i).chunks(80) {
                w.write_all(line)?;
                w.write_all(b"\n")?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_records() {
        let text = b"\n>chr1 first one\nACGT\nac\n>empty\n>chr2\r\nNNNN\r\nG\n";
        let mut r = FastaReader::new(&text[..]);
        let mut buf = Vec::new();
        let h = r.next_record(&mut buf).unwrap().unwrap();
        assert_eq!((h.name.as_str(), h.description.as_str(), &buf[..]), ("chr1", "first one", &b"ACGTac"[..]));
        let h = r.next_record(&mut buf).unwrap().unwrap();
        assert_eq!((h.name.as_str(), buf.len()), ("empty", 0));
        let h = r.next_record(&mut buf).unwrap().unwrap();
        assert_eq!((h.name.as_str(), &buf[..]), ("chr2", &b"NNNNG"[..]));
        assert!(r.next_record(&mut buf).unwrap().is_none());
    }

    #[test]
    fn genome_round_trip() {
        let mut g = Genome::default();
        g.push("a".into(), "desc".into(), &b"ACGT".repeat(50)).unwrap();
        g.push("b".into(), String::new(), b"NNacg").unwrap();
        let mut out = Vec::new();
        g.write_fasta(&mut out).unwrap();
        let mut r = FastaReader::new(&out[..]);
        let mut back = Genome::default();
        let mut buf = Vec::new();
        while let Some(h) = r.next_record(&mut buf).unwrap() {
            back.push(h.name, h.description, &buf).unwrap();
        }
        assert_eq!(back.contigs, g.contigs);
        assert_eq!(back.seq, g.seq);
    }
}
