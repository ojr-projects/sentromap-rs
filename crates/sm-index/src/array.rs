//! Flat typed arrays on disk: a 64-byte header followed by little-endian elements,
//! memory-mapped and used in place (design §6).

use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};
use std::marker::PhantomData;
use std::ops::Deref;
use std::path::Path;

use anyhow::{Context, Result, bail, ensure};
use bytemuck::{Pod, Zeroable};
use memmap2::{Mmap, MmapOptions};

pub const MAGIC: [u8; 8] = *b"SMAPARR\0";
pub const FORMAT_VERSION: u32 = 1;
pub const HEADER_LEN: usize = 64;

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable, Debug)]
struct Header {
    magic: [u8; 8],
    kind: [u8; 16],
    version: u32,
    elem_size: u32,
    count: u64,
    /// Kind-specific extra value (for bit vectors: the length in bits).
    aux: u64,
    _reserved: [u8; 16],
}
const _: () = assert!(size_of::<Header>() == HEADER_LEN);

fn kind_bytes(kind: &str) -> [u8; 16] {
    let mut k = [0u8; 16];
    assert!(kind.len() <= 16);
    k[..kind.len()].copy_from_slice(kind.as_bytes());
    k
}

/// Streams elements to disk, then patches the header with the final count.
pub struct ArrayWriter<T: Pod> {
    out: BufWriter<File>,
    kind: [u8; 16],
    count: u64,
    aux: u64,
    _t: PhantomData<T>,
}

impl<T: Pod> ArrayWriter<T> {
    pub fn create(path: &Path, kind: &str) -> Result<Self> {
        let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;
        let mut out = BufWriter::with_capacity(1 << 20, file);
        out.write_all(&[0u8; HEADER_LEN])?;
        Ok(Self { out, kind: kind_bytes(kind), count: 0, aux: 0, _t: PhantomData })
    }

    pub fn push(&mut self, x: T) -> Result<()> {
        self.out.write_all(bytemuck::bytes_of(&x))?;
        self.count += 1;
        Ok(())
    }

    pub fn extend(&mut self, xs: &[T]) -> Result<()> {
        self.out.write_all(bytemuck::cast_slice(xs))?;
        self.count += xs.len() as u64;
        Ok(())
    }

    pub fn len(&self) -> u64 {
        self.count
    }

    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    pub fn set_aux(&mut self, aux: u64) {
        self.aux = aux;
    }

    pub fn finish(mut self) -> Result<()> {
        let h = Header {
            magic: MAGIC,
            kind: self.kind,
            version: FORMAT_VERSION,
            elem_size: size_of::<T>() as u32,
            count: self.count,
            aux: self.aux,
            _reserved: [0; 16],
        };
        self.out.flush()?;
        let mut f = self.out.into_inner().map_err(|e| e.into_error())?;
        f.seek(SeekFrom::Start(0))?;
        f.write_all(bytemuck::bytes_of(&h))?;
        f.sync_all().ok();
        Ok(())
    }
}

/// Write a whole slice as an array file.
pub fn write_array<T: Pod>(path: &Path, kind: &str, xs: &[T], aux: u64) -> Result<()> {
    let mut w = ArrayWriter::create(path, kind)?;
    w.extend(xs)?;
    w.set_aux(aux);
    w.finish()
}

/// Options for mapping index files.
#[derive(Clone, Copy, Debug, Default)]
pub struct MapOptions {
    /// Pre-fault every page (`MAP_POPULATE`), for servers that want steady latency.
    pub populate: bool,
}

/// A memory-mapped array; dereferences to `&[T]`.
pub struct Array<T: Pod> {
    map: Mmap,
    len: usize,
    aux: u64,
    _t: PhantomData<T>,
}

impl<T: Pod> Array<T> {
    pub fn open(path: &Path, kind: &str, opts: MapOptions) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let mut mo = MmapOptions::new();
        if opts.populate {
            mo.populate();
        }
        // SAFETY: index files are written once by the builder and treated as read-only.
        let map = unsafe { mo.map(&file) }.with_context(|| format!("mapping {}", path.display()))?;
        ensure!(map.len() >= HEADER_LEN, "{}: too short for a header", path.display());
        let h: Header = bytemuck::pod_read_unaligned(&map[..HEADER_LEN]);
        if h.magic != MAGIC {
            bail!("{}: not a sentromap array file", path.display());
        }
        ensure!(
            h.version == FORMAT_VERSION,
            "{}: format version {} (expected {FORMAT_VERSION})",
            path.display(),
            h.version
        );
        ensure!(h.kind == kind_bytes(kind), "{}: expected a {kind:?} array", path.display());
        ensure!(
            h.elem_size as usize == size_of::<T>(),
            "{}: element size {} (expected {})",
            path.display(),
            h.elem_size,
            size_of::<T>()
        );
        let len = h.count as usize;
        ensure!(map.len() == HEADER_LEN + len * size_of::<T>(), "{}: truncated or oversized", path.display());
        Ok(Self { map, len, aux: h.aux, _t: PhantomData })
    }

    /// Element size recorded in a file's header, without mapping it as a particular type.
    pub fn peek_elem_size(path: &Path) -> Result<u32> {
        let mut buf = [0u8; HEADER_LEN];
        std::io::Read::read_exact(&mut File::open(path)?, &mut buf)?;
        let h: Header = bytemuck::pod_read_unaligned(&buf);
        ensure!(h.magic == MAGIC, "{}: not a sentromap array file", path.display());
        Ok(h.elem_size)
    }

    pub fn aux(&self) -> u64 {
        self.aux
    }

    /// Advise the kernel about access patterns (e.g. sequential for scans).
    pub fn advise(&self, advice: memmap2::Advice) {
        let _ = self.map.advise(advice);
    }
}

impl<T: Pod> Deref for Array<T> {
    type Target = [T];
    #[inline]
    fn deref(&self) -> &[T] {
        bytemuck::cast_slice(&self.map[HEADER_LEN..HEADER_LEN + self.len * size_of::<T>()])
    }
}
