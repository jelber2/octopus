// Converted from C++ to Rust

use crate::core::types::haplotype::Haplotype;
use crate::basics::aligned_read::AlignedRead;
use crate::basics::genomic_region::GenomicRegion;

pub struct BamRealigner {
    min_realignment_score: f64,
}

impl BamRealigner {
    pub fn new(min_realignment_score: f64) -> Self {
        BamRealigner { min_realignment_score }
    }

    pub fn realign(
        &self,
        reads: &mut Vec<AlignedRead>,
        haplotypes: &[Haplotype],
        region: &GenomicRegion,
    ) {
    }
}
