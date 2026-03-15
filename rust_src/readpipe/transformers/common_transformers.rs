// Converted from C++ to Rust

use crate::basics::aligned_read::AlignedRead;
use super::read_transformer::ReadTransformer;

pub struct CapBaseQualitiesTransformer {
    max_quality: u8,
}

impl CapBaseQualitiesTransformer {
    pub fn new(max_quality: u8) -> Self { CapBaseQualitiesTransformer { max_quality } }
}

impl ReadTransformer for CapBaseQualitiesTransformer {
    fn name(&self) -> &str { "CapBaseQualities" }

    fn transform(&self, read: &mut AlignedRead) {
        for q in read.qualities_mut() {
            if *q > self.max_quality {
                *q = self.max_quality;
            }
        }
    }
}
