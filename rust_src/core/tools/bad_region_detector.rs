// Converted from C++ to Rust

use crate::basics::genomic_region::GenomicRegion;
use crate::basics::aligned_read::AlignedRead;

pub struct BadRegionDetector {
    max_coverage: usize,
    min_coverage: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum RegionStatus {
    Good,
    HighCoverage,
    LowCoverage,
    HighMismatchRate,
}

impl BadRegionDetector {
    pub fn new(min_coverage: usize, max_coverage: usize) -> Self {
        BadRegionDetector { min_coverage, max_coverage }
    }

    pub fn is_bad(&self, reads: &[AlignedRead], region: &GenomicRegion) -> RegionStatus {
        let coverage = reads.len();

        if coverage > self.max_coverage {
            RegionStatus::HighCoverage
        } else if coverage < self.min_coverage {
            RegionStatus::LowCoverage
        } else {
            RegionStatus::Good
        }
    }
}
