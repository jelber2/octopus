// Converted from C++ to Rust

use crate::basics::aligned_read::AlignedRead;
use crate::basics::genomic_region::GenomicRegion;

#[derive(Debug, Clone, Default)]
pub struct ReadSetProfile {
    pub mean_read_length: f64,
    pub read_length_std: f64,
    pub mean_insert_size: f64,
    pub insert_size_std: f64,
    pub mean_coverage: f64,
}

pub fn profile_reads(reads: &[AlignedRead], _region: &GenomicRegion) -> ReadSetProfile {
    if reads.is_empty() {
        return ReadSetProfile::default();
    }

    let lengths: Vec<f64> = reads.iter().map(|r| r.sequence().len() as f64).collect();
    let mean_len = lengths.iter().sum::<f64>() / lengths.len() as f64;
    let var_len = lengths.iter().map(|&l| (l - mean_len).powi(2)).sum::<f64>() / lengths.len() as f64;

    ReadSetProfile {
        mean_read_length: mean_len,
        read_length_std: var_len.sqrt(),
        mean_insert_size: 0.0,
        insert_size_std: 0.0,
        mean_coverage: 0.0,
    }
}
