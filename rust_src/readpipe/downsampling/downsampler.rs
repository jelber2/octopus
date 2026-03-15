// Converted from C++ to Rust

use crate::basics::aligned_read::AlignedRead;
use crate::basics::genomic_region::GenomicRegion;

pub trait Downsampler: Send + Sync {
    fn downsample(&self, reads: &mut Vec<AlignedRead>, region: &GenomicRegion);
    fn name(&self) -> &str;
}
