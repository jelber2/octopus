// Converted from C++ to Rust

use crate::basics::aligned_read::AlignedRead;
use crate::basics::genomic_region::GenomicRegion;
use super::downsampler::Downsampler;

pub struct CoverageDownsampler {
    max_coverage: usize,
}

impl CoverageDownsampler {
    pub fn new(max_coverage: usize) -> Self {
        CoverageDownsampler { max_coverage }
    }
}

impl Downsampler for CoverageDownsampler {
    fn name(&self) -> &str { "CoverageDownsampler" }

    fn downsample(&self, reads: &mut Vec<AlignedRead>, region: &GenomicRegion) {
        if reads.len() <= self.max_coverage { return; }

        let step = reads.len() / self.max_coverage;
        let mut kept: Vec<AlignedRead> = reads.iter().enumerate()
            .filter(|(i, _)| i % step == 0)
            .take(self.max_coverage)
            .map(|(_, r)| r.clone())
            .collect();

        *reads = kept;
    }
}
