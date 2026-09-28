//! Core types for sentromap: k-mer encoding and distances, the global coordinate space,
//! FASTA input, the brute-force oracle and synthetic test genomes.

pub mod ball;
pub mod contig;
pub mod fasta;
pub mod kmer;
pub mod oracle;
pub mod synth;

pub use contig::{Contig, Contigs};
pub use kmer::{K, Kmer, canonical, diff_mask, hamming, revcomp, rotate, unrotate};
pub use oracle::{Hit, Strand};
