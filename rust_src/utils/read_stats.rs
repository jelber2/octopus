// Converted from C++ to Rust

use crate::basics::aligned_read::AlignedRead;

#[derive(Debug, Clone, Default)]
pub struct ReadStats {
    pub count: usize,
    pub mean_length: f64,
    pub mean_mapping_quality: f64,
    pub duplicate_fraction: f64,
}

pub fn compute_stats(reads: &[AlignedRead]) -> ReadStats {
    if reads.is_empty() {
        return ReadStats::default();
    }
    let count = reads.len();
    let total_length: usize = reads.iter().map(|r| r.sequence().len()).sum();
    let total_mq: f64 = reads.iter().map(|r| r.mapping_quality() as f64).sum();
    let duplicates = reads.iter().filter(|r| r.is_duplicate()).count();

    ReadStats {
        count,
        mean_length: total_length as f64 / count as f64,
        mean_mapping_quality: total_mq / count as f64,
        duplicate_fraction: duplicates as f64 / count as f64,
    }
}
